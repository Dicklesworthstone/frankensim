//! 3D Morton (Z-order) codec — the spatial identity underlying tile layout
//! (plan §5.3): Z-order tile ranks are simultaneously stencil-friendly,
//! sparse-friendly (FrankenVDB leaves align to them), and give spatial work
//! stable logical identities for scheduling, RNG keying, and reductions.
//!
//! Backends: portable magic-bits (the correctness reference, all targets)
//! and a BMI2 PDEP/PEXT capsule on x86-64 (`bmi2/`), selected ONCE into a
//! function table per the dispatch doctrine (plan §5.1 consequence 5).
//! G0 law: the backends are bit-identical over the whole domain.

use std::sync::OnceLock;

#[cfg(target_arch = "x86_64")]
mod bmi2;

/// Bits per axis: coordinates must be `< 2^21` so three interleaved axes fit
/// one `u64` Morton code.
pub const MORTON_BITS: u32 = 21;

/// Exclusive upper bound for each coordinate (`2^21`).
pub const MORTON_COORD_LIMIT: u32 = 1 << MORTON_BITS;

struct MortonFns {
    encode: fn(u32, u32, u32) -> u64,
    decode: fn(u64) -> (u32, u32, u32),
    backend: &'static str,
}

static FNS: OnceLock<MortonFns> = OnceLock::new();

fn fns() -> &'static MortonFns {
    FNS.get_or_init(|| {
        #[cfg(target_arch = "x86_64")]
        if std::arch::is_x86_feature_detected!("bmi2") {
            return MortonFns {
                encode: bmi2::encode,
                decode: bmi2::decode,
                backend: "bmi2",
            };
        }
        MortonFns {
            encode: encode_magic,
            decode: decode_magic,
            backend: "magic-bits",
        }
    })
}

/// The backend selected for this process ("bmi2" or "magic-bits") — a
/// ledger/log fact, resolved exactly once.
#[must_use]
pub fn morton_backend() -> &'static str {
    fns().backend
}

/// Interleave `(x, y, z)` into a Morton code. Coordinates above
/// [`MORTON_COORD_LIMIT`] are masked to their low 21 bits (callers —
/// `TileGrid` — enforce domain bounds with a structured error first).
#[inline]
#[must_use]
pub fn morton3_encode(x: u32, y: u32, z: u32) -> u64 {
    (fns().encode)(x, y, z)
}

/// Invert [`morton3_encode`]. Bits above `3 * MORTON_BITS` are ignored.
#[inline]
#[must_use]
pub fn morton3_decode(code: u64) -> (u32, u32, u32) {
    (fns().decode)(code)
}

/// Spread the low 21 bits of `v` so bit i lands at bit 3i (magic-bits
/// reference; constants are the canonical 64-bit 3D Morton set).
fn part1by2(v: u64) -> u64 {
    let mut v = v & 0x1f_ffff;
    v = (v | (v << 32)) & 0x1f_0000_0000_ffff;
    v = (v | (v << 16)) & 0x1f_0000_ff00_00ff;
    v = (v | (v << 8)) & 0x100f_00f0_0f00_f00f;
    v = (v | (v << 4)) & 0x10c3_0c30_c30c_30c3;
    v = (v | (v << 2)) & 0x1249_2492_4924_9249;
    v
}

/// Inverse of [`part1by2`]: gather every third bit back into the low 21.
fn compact1by2(v: u64) -> u64 {
    let mut v = v & 0x1249_2492_4924_9249;
    v = (v ^ (v >> 2)) & 0x10c3_0c30_c30c_30c3;
    v = (v ^ (v >> 4)) & 0x100f_00f0_0f00_f00f;
    v = (v ^ (v >> 8)) & 0x1f_0000_ff00_00ff;
    v = (v ^ (v >> 16)) & 0x1f_0000_0000_ffff;
    v = (v ^ (v >> 32)) & 0x1f_ffff;
    v
}

pub(crate) fn encode_magic(x: u32, y: u32, z: u32) -> u64 {
    part1by2(u64::from(x)) | (part1by2(u64::from(y)) << 1) | (part1by2(u64::from(z)) << 2)
}

pub(crate) fn decode_magic(code: u64) -> (u32, u32, u32) {
    (
        compact1by2(code) as u32,
        compact1by2(code >> 1) as u32,
        compact1by2(code >> 2) as u32,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    /// In-house LCG (same constants as the fs-qty hardening battery).
    struct Lcg(u64);
    impl Lcg {
        fn next(&mut self) -> u64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            self.0
        }
    }

    #[test]
    fn known_answers_anchor_the_bit_order() {
        // x occupies bit 0, y bit 1, z bit 2 (x is the fastest axis).
        assert_eq!(encode_magic(1, 0, 0), 0b001);
        assert_eq!(encode_magic(0, 1, 0), 0b010);
        assert_eq!(encode_magic(0, 0, 1), 0b100);
        assert_eq!(encode_magic(3, 0, 0), 0b001_001);
        assert_eq!(encode_magic(7, 7, 7), 0b111_111_111);
        assert_eq!(
            encode_magic(MORTON_COORD_LIMIT - 1, 0, 0),
            0x1249_2492_4924_9249
        );
    }

    #[test]
    fn magic_bits_bijection_exhaustive_small_and_random_large() {
        for x in 0..16 {
            for y in 0..16 {
                for z in 0..16 {
                    let c = encode_magic(x, y, z);
                    assert_eq!(decode_magic(c), (x, y, z));
                }
            }
        }
        let mut rng = Lcg(0x5EED_0001);
        for _ in 0..100_000 {
            let x = (rng.next() as u32) & (MORTON_COORD_LIMIT - 1);
            let y = (rng.next() as u32) & (MORTON_COORD_LIMIT - 1);
            let z = (rng.next() as u32) & (MORTON_COORD_LIMIT - 1);
            let c = encode_magic(x, y, z);
            assert_eq!(decode_magic(c), (x, y, z), "seed case ({x},{y},{z})");
            assert!(c < 1 << (3 * MORTON_BITS));
        }
    }

    #[test]
    fn dispatched_backend_matches_the_magic_reference() {
        // G0 tier-equivalence law (exercises the BMI2 capsule on x86 hosts
        // that have it; degenerate-but-valid on aarch64 where the dispatch
        // IS the reference).
        let mut rng = Lcg(0xB141_2026);
        for _ in 0..100_000 {
            let x = (rng.next() as u32) & (MORTON_COORD_LIMIT - 1);
            let y = (rng.next() as u32) & (MORTON_COORD_LIMIT - 1);
            let z = (rng.next() as u32) & (MORTON_COORD_LIMIT - 1);
            assert_eq!(morton3_encode(x, y, z), encode_magic(x, y, z));
            let c = encode_magic(x, y, z);
            assert_eq!(morton3_decode(c), decode_magic(c));
        }
        assert!(matches!(morton_backend(), "bmi2" | "magic-bits"));
    }

    #[test]
    fn zorder_is_monotone_within_an_octant_row() {
        // Locality sanity: consecutive x within an aligned 2-block differ
        // only in the low bit of the code.
        let a = encode_magic(4, 6, 2);
        let b = encode_magic(5, 6, 2);
        assert_eq!(a ^ b, 1);
    }
}
