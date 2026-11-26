//! Canonical replay state for the full-width executor stream.

use crate::wide::{STREAM_SEMANTICS_VERSION_V2, StreamKeyV2, StreamV2};
use crate::{STREAM_CHECKPOINT_MAGIC, StreamReplayError};

/// Version of the Philox4x64 checkpoint frame.
pub const STREAM_CHECKPOINT_VERSION_V2: u32 = 2;
/// Identity domain included verbatim in the v2 frame.
pub const STREAM_CHECKPOINT_IDENTITY_DOMAIN_V2: &str =
    "org.frankensim.fs-rand.stream-checkpoint.v2";

const MAGIC_END: usize = STREAM_CHECKPOINT_MAGIC.len();
const DOMAIN_END: usize = MAGIC_END + STREAM_CHECKPOINT_IDENTITY_DOMAIN_V2.len();
const VERSION_END: usize = DOMAIN_END + 4;
const SEMANTICS_END: usize = VERSION_END + 4;
const POSITION_END: usize = SEMANTICS_END + 5 * 8;
/// Fixed frame width including the trailing little-endian FNV-1a checksum.
pub const STREAM_CHECKPOINT_CANONICAL_LEN_V2: usize = POSITION_END + 8;

/// Full executor stream identity and next draw position for replay.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct StreamCheckpointV2 {
    /// Retained transport version.
    pub checkpoint_version: u32,
    /// Version of the mapping from a position to its random bits.
    pub stream_semantics_version: u32,
    /// All four u64 executor identity fields.
    pub key: StreamKeyV2,
    /// Next draw index.
    pub index: u64,
}

/// Replay admission failure for a full-width checkpoint.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StreamReplayErrorV2 {
    /// Invalid framing or an unsupported transport or draw version.
    Format(StreamReplayError),
    /// The stored FNV-1a checksum differs from the checksum of the frame.
    ChecksumMismatch,
}

impl core::fmt::Display for StreamReplayErrorV2 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Format(error) => core::fmt::Display::fmt(error, f),
            Self::ChecksumMismatch => f.write_str("stream checkpoint checksum mismatch"),
        }
    }
}

impl core::error::Error for StreamReplayErrorV2 {
    fn source(&self) -> Option<&(dyn core::error::Error + 'static)> {
        match self {
            Self::Format(error) => Some(error),
            Self::ChecksumMismatch => None,
        }
    }
}

impl From<StreamReplayError> for StreamReplayErrorV2 {
    fn from(error: StreamReplayError) -> Self {
        Self::Format(error)
    }
}

// FNV-1a over the preceding bytes detects accidental corruption. It plays no
// part in stream identity and does not authenticate untrusted input.
fn frame_checksum(bytes: &[u8]) -> u64 {
    bytes.iter().fold(0xcbf2_9ce4_8422_2325, |hash, &byte| {
        (hash ^ u64::from(byte)).wrapping_mul(0x0100_0000_01b3)
    })
}

impl StreamCheckpointV2 {
    /// Capture a position under the current full-width transport and semantics.
    #[must_use]
    pub const fn current(key: StreamKeyV2, index: u64) -> Self {
        Self {
            checkpoint_version: STREAM_CHECKPOINT_VERSION_V2,
            stream_semantics_version: STREAM_SEMANTICS_VERSION_V2,
            key,
            index,
        }
    }

    /// Encode magic, domain, both versions, seed, kernel, tile, iteration,
    /// next index, and a checksum of all preceding bytes, in that order.
    /// All integers use little-endian byte order.
    #[must_use]
    pub fn to_canonical_le_bytes(self) -> [u8; STREAM_CHECKPOINT_CANONICAL_LEN_V2] {
        let mut bytes = [0; STREAM_CHECKPOINT_CANONICAL_LEN_V2];
        bytes[..MAGIC_END].copy_from_slice(&STREAM_CHECKPOINT_MAGIC);
        bytes[MAGIC_END..DOMAIN_END]
            .copy_from_slice(STREAM_CHECKPOINT_IDENTITY_DOMAIN_V2.as_bytes());
        bytes[DOMAIN_END..VERSION_END].copy_from_slice(&self.checkpoint_version.to_le_bytes());
        bytes[VERSION_END..SEMANTICS_END]
            .copy_from_slice(&self.stream_semantics_version.to_le_bytes());
        let StreamKeyV2 {
            seed,
            kernel,
            tile,
            iteration,
        } = self.key;
        for (dest, value) in bytes[SEMANTICS_END..POSITION_END]
            .chunks_exact_mut(8)
            .zip([seed, kernel, tile, iteration, self.index])
        {
            dest.copy_from_slice(&value.to_le_bytes());
        }
        let checksum = frame_checksum(&bytes[..POSITION_END]);
        bytes[POSITION_END..].copy_from_slice(&checksum.to_le_bytes());
        bytes
    }

    /// Decode exactly one current frame, checking both versions and checksum.
    pub fn from_canonical_le_bytes(bytes: &[u8]) -> Result<Self, StreamReplayErrorV2> {
        if bytes.len() != STREAM_CHECKPOINT_CANONICAL_LEN_V2 {
            return Err(StreamReplayError::InvalidCheckpointLength {
                actual: bytes.len(),
                expected: STREAM_CHECKPOINT_CANONICAL_LEN_V2,
            }
            .into());
        }
        if bytes[..MAGIC_END] != STREAM_CHECKPOINT_MAGIC {
            return Err(StreamReplayError::InvalidCheckpointMagic.into());
        }
        if &bytes[MAGIC_END..DOMAIN_END] != STREAM_CHECKPOINT_IDENTITY_DOMAIN_V2.as_bytes() {
            return Err(StreamReplayError::InvalidCheckpointDomain.into());
        }
        let [seed, kernel, tile, iteration, index] = core::array::from_fn(|i| {
            let start = SEMANTICS_END + i * 8;
            u64::from_le_bytes(
                bytes[start..start + 8]
                    .try_into()
                    .expect("fixed-width field"),
            )
        });
        let checkpoint = Self {
            checkpoint_version: u32::from_le_bytes(
                bytes[DOMAIN_END..VERSION_END]
                    .try_into()
                    .expect("four-byte version"),
            ),
            stream_semantics_version: u32::from_le_bytes(
                bytes[VERSION_END..SEMANTICS_END]
                    .try_into()
                    .expect("four-byte version"),
            ),
            key: StreamKeyV2 {
                seed,
                kernel,
                tile,
                iteration,
            },
            index,
        };
        checkpoint.validate_versions()?;
        if bytes[POSITION_END..] != frame_checksum(&bytes[..POSITION_END]).to_le_bytes() {
            return Err(StreamReplayErrorV2::ChecksumMismatch);
        }
        Ok(checkpoint)
    }

    fn validate_versions(self) -> Result<(), StreamReplayError> {
        if self.checkpoint_version != STREAM_CHECKPOINT_VERSION_V2 {
            return Err(StreamReplayError::UnknownCheckpointVersion {
                declared: self.checkpoint_version,
                supported: STREAM_CHECKPOINT_VERSION_V2,
            });
        }
        if self.stream_semantics_version != STREAM_SEMANTICS_VERSION_V2 {
            return Err(StreamReplayError::UnknownStreamSemanticsVersion {
                declared: self.stream_semantics_version,
                supported: STREAM_SEMANTICS_VERSION_V2,
            });
        }
        Ok(())
    }
}

impl StreamV2 {
    /// Capture the full identity and next draw index.
    #[must_use]
    pub const fn checkpoint(&self) -> StreamCheckpointV2 {
        StreamCheckpointV2::current(self.key, self.index)
    }

    /// Resume a checkpoint carrying the exact v2 transport and draw versions.
    pub fn resume(checkpoint: StreamCheckpointV2) -> Result<Self, StreamReplayErrorV2> {
        checkpoint.validate_versions()?;
        Ok(Self {
            key: checkpoint.key,
            index: checkpoint.index,
        })
    }

    /// Decode and resume a canonical v2 checkpoint.
    pub fn resume_retained(bytes: &[u8]) -> Result<Self, StreamReplayErrorV2> {
        Self::resume(StreamCheckpointV2::from_canonical_le_bytes(bytes)?)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{StreamCheckpoint, StreamKey};

    const KEY: StreamKeyV2 = StreamKeyV2 {
        seed: 0x0807_0605_0403_0201,
        kernel: 0x100f_0e0d_0c0b_0a09,
        tile: 0x1817_1615_1413_1211,
        iteration: 0x201f_1e1d_1c1b_1a19,
    };
    const INDEX: u64 = 0x2827_2625_2423_2221;

    #[test]
    fn canonical_frame_and_continuation() {
        let checkpoint = StreamCheckpointV2::current(KEY, INDEX);
        let bytes = checkpoint.to_canonical_le_bytes();
        let expected = *b"FSRCKPT\0org.frankensim.fs-rand.stream-checkpoint.v2\
            \x02\0\0\0\x02\0\0\0\
            \x01\x02\x03\x04\x05\x06\x07\x08\x09\x0a\x0b\x0c\x0d\x0e\x0f\x10\
            \x11\x12\x13\x14\x15\x16\x17\x18\x19\x1a\x1b\x1c\x1d\x1e\x1f\x20\
            \x21\x22\x23\x24\x25\x26\x27\x28\x3e\x71\xbc\x62\x2c\x0e\x05\x29";
        assert_eq!(bytes, expected);
        assert_eq!(
            StreamCheckpointV2::from_canonical_le_bytes(&bytes),
            Ok(checkpoint)
        );
        let mut stream = StreamV2::resume_retained(&bytes).unwrap();
        assert_eq!(stream.checkpoint(), checkpoint);
        for offset in 0..8 {
            assert_eq!(stream.next_u64(), StreamV2::at(KEY, INDEX + offset)[0]);
        }
    }

    #[test]
    fn damaged_and_cross_version_frames_are_refused() {
        let bytes = StreamCheckpointV2::current(KEY, u64::MAX).to_canonical_le_bytes();
        for i in 0..bytes.len() {
            let mut damaged = bytes;
            damaged[i] ^= 1;
            assert!(StreamV2::resume_retained(&damaged).is_err(), "byte {i}");
        }
        assert!(StreamV2::resume_retained(&bytes[..bytes.len() - 1]).is_err());
        let mut extended = bytes.to_vec();
        extended.push(0);
        assert!(StreamV2::resume_retained(&extended).is_err());

        let legacy = StreamCheckpoint::current(
            StreamKey {
                seed: 1,
                kernel: 2,
                tile: 3,
            },
            4,
        );
        assert!(StreamV2::resume_retained(&legacy.to_canonical_le_bytes()).is_err());
        assert!(StreamCheckpoint::from_canonical_le_bytes(&bytes).is_err());
    }

    #[test]
    fn stale_or_unknown_versions_are_rejected_even_with_a_valid_checksum() {
        let current = StreamCheckpointV2::current(KEY, INDEX);
        for version in [1, 3, u32::MAX] {
            for checkpoint in [
                StreamCheckpointV2 {
                    checkpoint_version: version,
                    ..current
                },
                StreamCheckpointV2 {
                    stream_semantics_version: version,
                    ..current
                },
            ] {
                assert!(StreamV2::resume(checkpoint).is_err());
                assert!(StreamV2::resume_retained(&checkpoint.to_canonical_le_bytes()).is_err());
            }
        }
    }
}
