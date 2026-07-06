//! fs-simd — SIMD tiers behind safe façades (plan §5.1, patch Rev Q):
//! Tier 0 scalar stable Rust (the correctness reference, always available);
//! Tier 1 `std::arch` leaf capsules — NEON (aarch64) and AVX2/AVX-512
//! (x86-64) — each a registered unsafe capsule with a SAFETY.md;
//! Tier 2 nightly portable-SIMD, feature-gated, never load-bearing.
//!
//! Dispatch: resolved ONCE into a function table ([`ops`]), keyed by
//! fs-substrate's tier detection — no per-call branching in hot loops.
//! Under Miri the table routes to scalar (capsule intrinsics are outside
//! Miri's model; the SAFETY.md files document the compensating checks).
//!
//! Determinism contract: per tier, every primitive has a FIXED evaluation /
//! reduction shape (same input → same bits on the same tier). ACROSS tiers,
//! elementwise fused ops match bitwise (FMA policy: scalar twin uses
//! `mul_add`); reductions may differ within a documented envelope — that
//! difference is machine identity (G5's cross-ISA report), never run jitter.

pub mod scalar;

#[cfg(all(target_arch = "aarch64", not(miri)))]
pub mod neon;

#[cfg(all(target_arch = "x86_64", not(miri)))]
pub mod x86;

use fs_substrate::SimdTier;
use std::sync::OnceLock;

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Ternary elementwise kernel signature (a, b, c, out).
pub type TernaryOp = fn(&[f64], &[f64], &[f64], &mut [f64]);

/// The resolved-once function table (plan §5.1 consequence 5).
pub struct Ops {
    /// Tier the table was built for (ledger/tune-table key material).
    pub tier: SimdTier,
    /// y[i] = a·x[i] + y[i] (fused).
    pub axpy: fn(f64, &[f64], &mut [f64]),
    /// x[i] *= a.
    pub scale: fn(f64, &mut [f64]),
    /// out[i] = a[i]·b[i].
    pub mul_elem: fn(&[f64], &[f64], &mut [f64]),
    /// out[i] = a[i]·b[i] + c[i] (fused).
    pub fma3: TernaryOp,
    /// Σ x[i]·y[i] (fixed per-tier shape).
    pub dot: fn(&[f64], &[f64]) -> f64,
    /// Σ x[i] (fixed per-tier shape).
    pub sum: fn(&[f64]) -> f64,
}

static OPS: OnceLock<Ops> = OnceLock::new();

/// The process-wide primitive table, resolved exactly once.
pub fn ops() -> &'static Ops {
    OPS.get_or_init(build_table)
}

const SCALAR_OPS: Ops = Ops {
    tier: SimdTier::Scalar,
    axpy: scalar::axpy,
    scale: scalar::scale,
    mul_elem: scalar::mul_elem,
    fma3: scalar::fma3,
    dot: scalar::dot,
    sum: scalar::sum,
};

fn build_table() -> Ops {
    #[cfg(miri)]
    {
        SCALAR_OPS
    }
    #[cfg(not(miri))]
    {
        match fs_substrate::dispatch_tier() {
            #[cfg(target_arch = "aarch64")]
            SimdTier::Neon => Ops {
                tier: SimdTier::Neon,
                axpy: neon::axpy,
                scale: neon::scale,
                mul_elem: neon::mul_elem,
                fma3: neon::fma3,
                dot: neon::dot,
                sum: neon::sum,
            },
            // x86 capsule v1 covers axpy/dot/sum (the <300-line capsule cap
            // is a feature: scale/mul_elem/fma3 arrive with their consumer,
            // fs-la's packing kernels). Fallbacks are the scalar twin.
            #[cfg(target_arch = "x86_64")]
            SimdTier::Avx2 | SimdTier::Avx512 => Ops {
                tier: fs_substrate::dispatch_tier(),
                axpy: x86::axpy,
                scale: scalar::scale,
                mul_elem: scalar::mul_elem,
                fma3: scalar::fma3,
                dot: x86::dot,
                sum: x86::sum,
            },
            _ => SCALAR_OPS,
        }
    }
}

/// True if `ptr` is aligned to the target's cache line (fs-substrate's
/// `CACHE_LINE`) — the padding/false-sharing audit helper.
#[must_use]
pub fn is_cache_line_aligned<T>(ptr: *const T) -> bool {
    (ptr as usize).is_multiple_of(fs_substrate::CACHE_LINE)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Deterministic input generator (LCG; fs-rand lands later in the graph).
    fn gen_vals(len: usize, seed: u64) -> Vec<f64> {
        let mut s = seed | 1;
        (0..len)
            .map(|i| {
                s = s
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                match (s >> 60) & 0x7 {
                    0 => 0.0,
                    1 => -0.0,
                    2 => f64::MIN_POSITIVE / 2.0, // subnormal
                    3 => 1e18, // large but products stay finite (envelope math needs finite)
                    _ => (((s >> 11) as f64) / (1u64 << 53) as f64 - 0.5) * (i as f64 + 1.0),
                }
            })
            .collect()
    }

    /// The battery both capsules cite in their SAFETY.md: every tail length,
    /// special values, elementwise-bitwise + reduction-envelope equivalence
    /// between the ACTIVE tier and the scalar twin.
    #[test]
    fn tier_equivalence_battery() {
        let t = ops();
        for len in 0..67 {
            for seed in [1u64, 42, 0xDEAD] {
                let x = gen_vals(len, seed);
                let y0 = gen_vals(len, seed ^ 0x7);
                let c = gen_vals(len, seed ^ 0x63);
                // axpy: bitwise (fused both sides).
                let mut y_tier = y0.clone();
                (t.axpy)(1.5, &x, &mut y_tier);
                let mut y_ref = y0.clone();
                scalar::axpy(1.5, &x, &mut y_ref);
                assert!(
                    y_tier
                        .iter()
                        .zip(&y_ref)
                        .all(|(a, b)| a.to_bits() == b.to_bits()),
                    "axpy diverged from twin at len {len} seed {seed} (tier {:?})",
                    t.tier
                );
                // scale: bitwise.
                let mut s_tier = x.clone();
                (t.scale)(-0.25, &mut s_tier);
                let mut s_ref = x.clone();
                scalar::scale(-0.25, &mut s_ref);
                assert!(
                    s_tier
                        .iter()
                        .zip(&s_ref)
                        .all(|(a, b)| a.to_bits() == b.to_bits())
                );
                // mul_elem / fma3: bitwise.
                let mut m_tier = vec![0.0; len];
                (t.mul_elem)(&x, &y0, &mut m_tier);
                let mut m_ref = vec![0.0; len];
                scalar::mul_elem(&x, &y0, &mut m_ref);
                assert!(
                    m_tier
                        .iter()
                        .zip(&m_ref)
                        .all(|(a, b)| a.to_bits() == b.to_bits())
                );
                let mut f_tier = vec![0.0; len];
                (t.fma3)(&x, &y0, &c, &mut f_tier);
                let mut f_ref = vec![0.0; len];
                scalar::fma3(&x, &y0, &c, &mut f_ref);
                assert!(
                    f_tier
                        .iter()
                        .zip(&f_ref)
                        .all(|(a, b)| a.to_bits() == b.to_bits())
                );
                // dot/sum: same-tier bit-stability + cross-shape envelope.
                let d1 = (t.dot)(&x, &y0);
                let d2 = (t.dot)(&x, &y0);
                assert_eq!(d1.to_bits(), d2.to_bits(), "same tier must be bit-stable");
                let d_ref = scalar::dot(&x, &y0);
                let scale_mag: f64 = x
                    .iter()
                    .zip(&y0)
                    .map(|(a, b)| (a * b).abs())
                    .sum::<f64>()
                    .max(1e-300);
                assert!(
                    (d1 - d_ref).abs() <= 1e-12 * scale_mag,
                    "dot outside envelope at len {len}: tier {d1} vs twin {d_ref}"
                );
                let s1 = (t.sum)(&x);
                let s_refv = scalar::sum(&x);
                let mag: f64 = x.iter().map(|v| v.abs()).sum::<f64>().max(1e-300);
                assert!((s1 - s_refv).abs() <= 1e-12 * mag);
            }
        }
        println!(
            "{{\"suite\":\"fs-simd/equivalence\",\"case\":\"battery\",\"verdict\":\"pass\",\"detail\":\"tier={} lens=0..67\"}}",
            t.tier.name()
        );
    }

    #[test]
    fn dispatch_table_is_singleton_and_tier_matches_substrate() {
        let a = std::ptr::from_ref(ops());
        let b = std::ptr::from_ref(ops());
        assert_eq!(a, b, "table must resolve once");
        #[cfg(all(target_arch = "aarch64", not(miri)))]
        assert_eq!(ops().tier, SimdTier::Neon);
        #[cfg(miri)]
        assert_eq!(ops().tier, SimdTier::Scalar);
    }

    #[test]
    fn known_answers_anchor_the_semantics() {
        // Small exact cases catch sign/lane-order bugs that equivalence
        // against a buggy twin could miss.
        let x = [1.0, 2.0, 3.0, 4.0, 5.0];
        let y = [10.0, 20.0, 30.0, 40.0, 50.0];
        assert_eq!((ops().dot)(&x, &y).to_bits(), 550.0f64.to_bits());
        assert_eq!((ops().sum)(&x).to_bits(), 15.0f64.to_bits());
        let mut z = y;
        (ops().axpy)(2.0, &x, &mut z);
        let want = [12.0f64, 24.0, 36.0, 48.0, 60.0];
        assert!(
            z.iter().zip(&want).all(|(a, b)| a.to_bits() == b.to_bits()),
            "{z:?}"
        );
    }

    #[test]
    fn cache_line_alignment_helper() {
        let v = vec![0u8; 256];
        let base = v.as_ptr() as usize;
        let aligned = base.next_multiple_of(fs_substrate::CACHE_LINE);
        assert!(is_cache_line_aligned(aligned as *const u8));
        assert!(!is_cache_line_aligned((aligned + 8) as *const u8));
    }

    #[test]
    #[should_panic(expected = "length mismatch")]
    fn length_mismatch_is_a_loud_programmer_error() {
        let x = [1.0, 2.0];
        let mut y = [1.0];
        (ops().axpy)(1.0, &x, &mut y);
    }
}
