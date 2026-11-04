//! Chunked acquisition with the original joint Welch segment owner.

use super::*;

/// Synchronized Welch estimation without retaining the full input record.
///
/// Declare the final sample count per channel at construction. This admits all
/// transforms and retained factors up front and gives every segment exactly
/// the same final-record normalization as [`WelchCrossSpectrum::estimate`].
/// Chunks are consecutive samples on the declared clock; channel identity and
/// order must remain fixed. Chunk boundaries never reset a window, overlap,
/// detrending or phase. There is no resampling, missing-data fill or padding.
///
/// Input storage is two segment-sized interleaved buffers: the overlapping
/// tail and its transactional replacement. Caller-owned chunks are borrowed;
/// completed input is discarded. The existing O(bins * segments * channels)
/// factor output, window, FFT buffers/plan and both input buffers are admitted
/// by `maximum_workspace_bytes`. No second factor bank or final transpose is
/// needed. Overlapping factors still imply neither independent realizations
/// nor cross-frequency coherence.
///
/// Every pushed chunk is atomic. Shape, length, nonfinite input, numerical
/// refusal or cancellation leaves accepted samples/segments unchanged, so the
/// caller may retry. Scratch and unaccepted future factor slots may change;
/// they are never exposed and are completely overwritten before acceptance.
pub struct WelchStream {
    output: WelchCrossSpectrum,
    transform: SegmentTransform,
    accepted_samples: usize,
    completed_segments: usize,
    pending: Vec<f64>,
    tail: Vec<f64>,
}

impl WelchStream {
    pub fn new(
        channels: usize,
        expected_samples_per_channel: usize,
        timestep_s: f64,
        config: WelchConfig,
        budget: WelchBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, WelchError> {
        poll(&mut cancelled)?;
        let shape = admit(
            channels,
            expected_samples_per_channel,
            timestep_s,
            config,
            budget,
            2,
        )?;
        let transform = SegmentTransform::new(timestep_s, config, shape.segments, &mut cancelled)?;
        let factors = allocate(shape.factor_count, C64::default())?;
        let mut pending = Vec::new();
        let mut tail = Vec::new();
        for buffer in [&mut pending, &mut tail] {
            buffer
                .try_reserve_exact(channels * config.segment_samples)
                .map_err(|_| WelchError::Allocation)?;
        }
        poll(&mut cancelled)?;
        Ok(Self {
            output: WelchCrossSpectrum {
                config,
                timestep_s,
                bin_width_hz: shape.bin_width_hz,
                channels,
                segments: shape.segments,
                input_samples: expected_samples_per_channel,
                used_samples: (shape.segments - 1) * config.hop_samples + config.segment_samples,
                factors,
            },
            transform,
            accepted_samples: 0,
            completed_segments: 0,
            pending,
            tail,
        })
    }

    /// Number of accepted samples in EACH channel, including an incomplete tail.
    pub const fn accepted_sample_count(&self) -> usize {
        self.accepted_samples
    }

    pub const fn remaining_sample_count(&self) -> usize {
        self.output.input_samples - self.accepted_samples
    }

    /// Complete windows belonging to accepted chunks, regardless of overlap.
    pub const fn completed_segment_count(&self) -> usize {
        self.completed_segments
    }

    /// Push equally long channel slices in their permanent physical order.
    /// An empty chunk is allowed. Even samples in the eventual discarded tail
    /// must be finite. A rejected chunk consumes no samples or segment slots.
    pub fn push(
        &mut self,
        channels: &[&[f64]],
        cancelled: impl FnMut() -> bool,
    ) -> Result<(), WelchError> {
        let samples = channels.first().map_or(0, |c| c.len());
        if channels.iter().any(|c| c.len() != samples) {
            return Err(WelchError::Invalid("channel sample counts differ"));
        }
        self.push_with(channels.len(), samples, |c, j| channels[c][j], cancelled)
    }

    /// Push synchronized Cartesian frames without a transpose or record copy.
    pub fn push_interleaved<const D: usize>(
        &mut self,
        samples: &[[f64; D]],
        cancelled: impl FnMut() -> bool,
    ) -> Result<(), WelchError> {
        self.push_with(D, samples.len(), |c, j| samples[j][c], cancelled)
    }

    fn push_with(
        &mut self,
        channels: usize,
        samples: usize,
        sample: impl Fn(usize, usize) -> f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), WelchError> {
        poll(&mut cancelled)?;
        if channels != self.output.channels || samples > self.remaining_sample_count() {
            return Err(WelchError::Invalid(
                "chunk exceeds the declared channel shape or record",
            ));
        }
        validate_samples(channels, samples, &sample, &mut cancelled)?;
        let accepted_samples = self.accepted_samples + samples;
        let config = self.output.config;
        let completed_segments = if accepted_samples < config.segment_samples {
            0
        } else {
            1 + (accepted_samples - config.segment_samples) / config.hop_samples
        };
        if completed_segments == self.completed_segments {
            // Small chunks append in place, avoiding a whole-window copy per
            // sample. No previously accepted element is overwritten.
            let old_len = self.pending.len();
            let append = (|| {
                for j in 0..samples {
                    for c in 0..channels {
                        if c % 256 == 0 {
                            poll(&mut cancelled)?;
                        }
                        self.pending.push(sample(c, j));
                    }
                }
                poll(&mut cancelled)
            })();
            if let Err(error) = append {
                self.pending.truncate(old_len);
                return Err(error);
            }
        } else {
            // Read a virtual concatenation of the old overlap and this borrowed
            // chunk. Large pushes need no input-sized allocation or copying.
            let pending_frames = self.pending.len() / channels;
            let read = |c, j| {
                if j < pending_frames {
                    self.pending[j * channels + c]
                } else {
                    sample(c, j - pending_frames)
                }
            };
            for segment in self.completed_segments..completed_segments {
                let start = (segment - self.completed_segments) * config.hop_samples;
                for channel in 0..channels {
                    self.transform.channel(
                        |j| read(channel, start + j),
                        |bin, value| {
                            self.output.factors
                                [(bin * self.output.segments + segment) * channels + channel] =
                                value;
                        },
                        &mut cancelled,
                    )?;
                }
            }
            self.tail.clear();
            let start = (completed_segments - self.completed_segments) * config.hop_samples;
            for j in start..pending_frames + samples {
                for c in 0..channels {
                    if c % 256 == 0 {
                        poll(&mut cancelled)?;
                    }
                    self.tail.push(read(c, j));
                }
            }
            // Commit only after every segment and the replacement tail exists.
            poll(&mut cancelled)?;
            core::mem::swap(&mut self.pending, &mut self.tail);
            self.tail.clear();
        }
        self.accepted_samples = accepted_samples;
        self.completed_segments = completed_segments;
        Ok(())
    }

    /// Consume the completed declared record and publish the original immutable
    /// spectrum, with its exact used-sample and discarded-tail counts. This is
    /// constant work: factors are already normalized and in their final layout.
    /// Check [`Self::remaining_sample_count`] before consuming an incomplete
    /// record; missing samples are refused, never padded or extrapolated.
    pub fn finish(self) -> Result<WelchCrossSpectrum, WelchError> {
        if self.remaining_sample_count() != 0 {
            Err(WelchError::Invalid(
                "declared synchronized record is incomplete",
            ))
        } else {
            Ok(self.output)
        }
    }
}
