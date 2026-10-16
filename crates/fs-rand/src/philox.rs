//! Philox4x32-10 (Random123 lineage): a COUNTER-BASED generator — the output
//! is a pure function of (key, counter), so parallel streams keyed by
//! LOGICAL work identity are reproducible regardless of which worker runs
//! which tile (Decalogue P2's RNG pillar, plan §6.7). Pure integer
//! arithmetic: bit-identical on every target by construction.

/// Philox 32-bit multipliers (Random123).
const M0: u32 = 0xD251_1F53;
const M1: u32 = 0xCD9E_8D57;
/// Weyl key-schedule increments (golden-ratio / sqrt(3)−1 fractions).
const W0: u32 = 0x9E37_79B9;
const W1: u32 = 0xBB67_AE85;

#[inline]
fn mulhilo(a: u32, b: u32) -> (u32, u32) {
    let wide = u64::from(a) * u64::from(b);
    ((wide >> 32) as u32, wide as u32)
}

/// One Philox round.
#[inline]
fn round(ctr: [u32; 4], key: [u32; 2]) -> [u32; 4] {
    let (hi0, lo0) = mulhilo(M0, ctr[0]);
    let (hi1, lo1) = mulhilo(M1, ctr[2]);
    [hi1 ^ ctr[1] ^ key[0], lo1, hi0 ^ ctr[3] ^ key[1], lo0]
}

/// The full 10-round Philox4x32 block function: 128-bit counter + 64-bit key
/// → 128 output bits.
#[must_use]
pub fn philox4x32_10(mut ctr: [u32; 4], mut key: [u32; 2]) -> [u32; 4] {
    for i in 0..10 {
        if i > 0 {
            key[0] = key[0].wrapping_add(W0);
            key[1] = key[1].wrapping_add(W1);
        }
        ctr = round(ctr, key);
    }
    ctr
}

/// BATCHED Philox4x32-10: `L` counter blocks under ONE shared key, computed
/// LANE-WISE (structure-of-arrays: `c[word][lane]`) so the integer rounds
/// auto-vectorize (SLP over the `L` lanes) for fill-heavy consumers (bead 1za9).
/// Each lane is BITWISE-IDENTICAL to the scalar [`philox4x32_10`] on the same
/// input — the fs-simd twin doctrine, verified in `tests/bulk.rs`. Pure integer
/// arithmetic, so the identity holds on every ISA regardless of whether the
/// compiler vectorized the lane loop.
///
/// A true hand-written NEON/AVX kernel (and its throughput gate) is a follow-up:
/// `fs-simd` exposes only float ops today and the workspace is on stable Rust
/// (no `portable_simd`); this SoA form is the correct, ISA-independent substrate
/// for that capsule.
#[must_use]
// The lane loops index four SoA arrays `c[word][lane]` by a shared `lane` — the
// vectorizable kernel form, not iterable as a single sequence.
#[allow(clippy::needless_range_loop)]
pub fn philox4x32_10_batch<const L: usize>(ctr: &[[u32; 4]; L], key: [u32; 2]) -> [[u32; 4]; L] {
    // Structure-of-arrays: c[word][lane].
    let mut c: [[u32; L]; 4] = [[0; L]; 4];
    for (lane, blk) in ctr.iter().enumerate() {
        for w in 0..4 {
            c[w][lane] = blk[w];
        }
    }
    // One shared running key (the fill case keys every block identically).
    let (mut k0, mut k1) = (key[0], key[1]);
    for r in 0..10 {
        if r > 0 {
            k0 = k0.wrapping_add(W0);
            k1 = k1.wrapping_add(W1);
        }
        for lane in 0..L {
            // All reads use the OLD lane state (computed before any write).
            let (hi0, lo0) = mulhilo(M0, c[0][lane]);
            let (hi1, lo1) = mulhilo(M1, c[2][lane]);
            let n0 = hi1 ^ c[1][lane] ^ k0;
            let n2 = hi0 ^ c[3][lane] ^ k1;
            c[0][lane] = n0;
            c[1][lane] = lo1;
            c[2][lane] = n2;
            c[3][lane] = lo0;
        }
    }
    let mut out = [[0u32; 4]; L];
    for (lane, o) in out.iter_mut().enumerate() {
        for w in 0..4 {
            o[w] = c[w][lane];
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The batched block function is bitwise-identical to the scalar per lane.
    #[test]
    fn batch_matches_scalar() {
        let key = [0xA409_3822, 0x299F_31D0];
        let ctr: [[u32; 4]; 8] =
            core::array::from_fn(|l| [0x243F_6A88 ^ l as u32, 0x85A3_08D3, l as u32, 7]);
        let batched = philox4x32_10_batch::<8>(&ctr, key);
        for (l, blk) in ctr.iter().enumerate() {
            assert_eq!(batched[l], philox4x32_10(*blk, key), "lane {l}");
        }
    }

    /// Known-answer tests from the Random123 distribution's kat_vectors
    /// (philox4x32, 10 rounds). These pin the EXACT semantics: any deviation
    /// in round structure, key schedule, or word order fails here.
    #[test]
    fn random123_known_answers() {
        let cases: [([u32; 4], [u32; 2], [u32; 4]); 3] = [
            (
                [0, 0, 0, 0],
                [0, 0],
                [0x6627_e8d5, 0xe169_c58d, 0xbc57_ac4c, 0x9b00_dbd8],
            ),
            (
                [0xffff_ffff; 4],
                [0xffff_ffff; 2],
                [0x408f_276d, 0x41c8_3b0e, 0xa20b_c7c6, 0x6d54_51fd],
            ),
            (
                // π digits as counter and key (the classic Random123 case).
                [0x243f_6a88, 0x85a3_08d3, 0x1319_8a2e, 0x0370_7344],
                [0xa409_3822, 0x299f_31d0],
                [0xd16c_fe09, 0x94fd_cceb, 0x5001_e420, 0x2412_6ea1],
            ),
        ];
        for (i, (ctr, key, want)) in cases.iter().enumerate() {
            let got = philox4x32_10(*ctr, *key);
            assert_eq!(
                got, *want,
                "KAT {i}: philox4x32-10({ctr:08x?}, {key:08x?}) = {got:08x?}, want {want:08x?}"
            );
        }
        println!(
            "{{\"suite\":\"fs-rand\",\"case\":\"philox-kat\",\"verdict\":\"pass\",\"detail\":\"3 Random123 vectors\"}}"
        );
    }

    #[test]
    fn counter_sensitivity_avalanche() {
        // Flipping any single counter bit must change roughly half the
        // output bits (crude avalanche check; the full statistical battery
        // is nightly-CI scope).
        let base = philox4x32_10([1, 2, 3, 4], [5, 6]);
        let mut total = 0u32;
        for word in 0..4 {
            for bit in 0..32 {
                let mut ctr = [1u32, 2, 3, 4];
                ctr[word] ^= 1 << bit;
                let out = philox4x32_10(ctr, [5, 6]);
                let diff: u32 = (0..4).map(|k| (out[k] ^ base[k]).count_ones()).sum();
                assert!(
                    (32..=96).contains(&diff),
                    "weak avalanche: ctr word {word} bit {bit} flipped only {diff}/128"
                );
                total += diff;
            }
        }
        let mean = f64::from(total) / 128.0;
        assert!(
            (54.0..=74.0).contains(&mean),
            "avalanche mean {mean} off 64"
        );
    }
}
