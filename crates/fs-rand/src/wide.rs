//! Full-width executor stream identities and Philox4x64-10 draw semantics.

/// Lossless bridge for the executor's four `u64` key fields.
pub const EXEC_KEY_BRIDGE_VERSION_V2: u32 = 2;

/// Philox4x64-10 with counter `[index, tile, iteration, kernel]`, key
/// `[seed, 2]`, and lane zero as the sequential draw.
pub const STREAM_SEMANTICS_VERSION_V2: u32 = 2;

/// Logical stream identity copied from an executor. Every field occupies
/// its own 64-bit word in the Philox position.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct StreamKeyV2 {
    /// Run seed.
    pub seed: u64,
    /// Executor kernel identifier.
    pub kernel: u64,
    /// Logical tile ordinal.
    pub tile: u64,
    /// Declared generation or iteration.
    pub iteration: u64,
}

impl StreamKeyV2 {
    /// Bridge the fields of `fs_exec::StreamKey` without narrowing or mixing.
    #[must_use]
    pub const fn from_exec_parts_v2(seed: u64, kernel_id: u64, tile: u64, iteration: u64) -> Self {
        Self {
            seed,
            kernel: kernel_id,
            tile,
            iteration,
        }
    }

    /// Start at draw index zero.
    #[must_use]
    pub const fn stream(self) -> StreamV2 {
        StreamV2 {
            key: self,
            index: 0,
        }
    }
}

/// A sequential cursor over Philox4x64-10 blocks.
#[derive(Debug, Clone, Copy)]
pub struct StreamV2 {
    pub(crate) key: StreamKeyV2,
    pub(crate) index: u64,
}

impl StreamV2 {
    /// Read all four output words at an absolute draw index.
    #[must_use]
    pub fn at(key: StreamKeyV2, index: u64) -> [u64; 4] {
        let StreamKeyV2 {
            seed,
            kernel,
            tile,
            iteration,
        } = key;
        philox4x64_10(
            [index, tile, iteration, kernel],
            [seed, u64::from(STREAM_SEMANTICS_VERSION_V2)],
        )
    }

    /// Index of the next sequential draw.
    #[must_use]
    pub const fn index(&self) -> u64 {
        self.index
    }

    /// Return word zero at the current index and advance modulo 2^64.
    #[must_use]
    pub fn next_u64(&mut self) -> u64 {
        let word = Self::at(self.key, self.index)[0];
        self.index = self.index.wrapping_add(1);
        word
    }

    /// Uniform in [0, 1) from the high 53 bits of one draw.
    #[must_use]
    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / 9_007_199_254_740_992.0)
    }
}

/// Ten rounds of Philox4x64 using the Random123 multipliers and Weyl increments.
fn philox4x64_10(mut counter: [u64; 4], mut key: [u64; 2]) -> [u64; 4] {
    const M0: u64 = 0xd2e7_470e_e14c_6c93;
    const M1: u64 = 0xca5a_8263_9512_1157;
    const W0: u64 = 0x9e37_79b9_7f4a_7c15;
    const W1: u64 = 0xbb67_ae85_84ca_a73b;

    for _ in 0..10 {
        let p0 = u128::from(M0) * u128::from(counter[0]);
        let p1 = u128::from(M1) * u128::from(counter[2]);
        counter = [
            (p1 >> 64) as u64 ^ counter[1] ^ key[0],
            p1 as u64,
            (p0 >> 64) as u64 ^ counter[3] ^ key[1],
            p0 as u64,
        ];
        key[0] = key[0].wrapping_add(W0);
        key[1] = key[1].wrapping_add(W1);
    }
    counter
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    #[test]
    fn random123_known_answer() {
        // DEShawResearch/random123, tests/kat_vectors, philox4x64 10 (pi digits).
        assert_eq!(
            philox4x64_10(
                [
                    0x243f_6a88_85a3_08d3,
                    0x1319_8a2e_0370_7344,
                    0xa409_3822_299f_31d0,
                    0x082e_fa98_ec4e_6c89
                ],
                [0x4528_21e6_38d0_1377, 0xbe54_66cf_34e9_0c6c],
            ),
            [
                0xa528_f454_03e6_1d95,
                0x38c7_2dbd_566e_9788,
                0xa5a1_610e_72fd_18b5,
                0x57bd_43b5_e52b_7fe6
            ],
        );
    }

    const WIDE_STREAM_V2_100K_FNV1A64_GOLDEN: u64 = 0x17b7_7422_223d_6ca0;

    #[test]
    fn wide_stream_golden() {
        // Independently generated with NumPy's Philox4x64, consuming lane
        // zero of each block and hashing its eight little-endian bytes.
        let key = StreamKeyV2::from_exec_parts_v2(
            0x0123_4567_89ab_cdef,
            0xfedc_ba98_7654_3210,
            (1 << 32) | 5,
            u64::MAX,
        );
        let mut stream = key.stream();
        let mut digest = 0xcbf2_9ce4_8422_2325_u64;
        for _ in 0..100_000 {
            for byte in stream.next_u64().to_le_bytes() {
                digest = (digest ^ u64::from(byte)).wrapping_mul(0x0100_0000_01b3);
            }
        }
        assert_eq!(digest, WIDE_STREAM_V2_100K_FNV1A64_GOLDEN);
        assert_eq!(stream.index(), 100_000);
    }

    #[test]
    fn executor_boundaries_preserve_all_four_fields() {
        let boundary = [
            0,
            1,
            5,
            u64::from(u32::MAX),
            1 << 32,
            (1 << 32) | 5,
            u64::MAX,
        ];
        let mut keys = HashSet::new();
        let mut first_draws = HashSet::new();
        for seed in boundary {
            for kernel in boundary {
                for tile in boundary {
                    for iteration in boundary {
                        let key = StreamKeyV2::from_exec_parts_v2(seed, kernel, tile, iteration);
                        assert_eq!(
                            [key.seed, key.kernel, key.tile, key.iteration],
                            [seed, kernel, tile, iteration]
                        );
                        assert!(keys.insert(key));
                        assert!(first_draws.insert(key.stream().next_u64()), "{key:?}");
                    }
                }
            }
        }
    }

    #[test]
    fn cursor_reads_absolute_indices_through_wrap() {
        let key = StreamKeyV2::from_exec_parts_v2(1, 2, 3, 4);
        let mut stream = StreamV2 {
            key,
            index: u64::MAX,
        };
        // NumPy Philox at counter [u64::MAX, 3, 4, 2] with key [1, 2].
        let block = [
            0x5686_d496_1f66_b652,
            0x1dfd_a065_4590_7457,
            0xf8b1_e599_258e_09b2,
            0x5a38_9fea_ba74_f059,
        ];
        assert_eq!(StreamV2::at(key, u64::MAX), block);
        assert_eq!(stream.next_u64(), block[0]);
        assert_eq!(stream.index(), 0);
        assert_eq!(stream.next_u64(), StreamV2::at(key, 0)[0]);
    }
}
