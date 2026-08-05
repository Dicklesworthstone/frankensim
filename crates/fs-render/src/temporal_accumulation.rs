//! Deterministic transactional accumulation across a resolved shutter.
//!
//! Samples are evaluated strictly in ascending absolute logical-sample order.
//! Each requested contiguous range is staged privately and committed only after
//! every sample and the final cancellation checkpoint succeed.  Splitting one
//! range into contiguous progressive partitions therefore executes the exact
//! same floating-point additions in the exact same order.

use core::{fmt, ops::Range};

use fs_exec::{Cancelled, Cx};

use crate::motion::{NormalizedShutterTime, ShutterInterval};

/// Interpretation of the accumulator's three linear channels.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TemporalColorSpace {
    /// Linear red, green, and blue channels.
    LinearRgb,
    /// CIE-style X, Y, and Z tristimulus channels.
    Xyz,
}

/// One deterministic evaluator request inside an admitted shutter.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TemporalSamplePoint {
    logical_sample_id: u64,
    normalized_time: NormalizedShutterTime,
    absolute_time_s: f64,
}

impl TemporalSamplePoint {
    /// Absolute logical sample identity used by deterministic shutter sampling.
    #[must_use]
    pub const fn logical_sample_id(self) -> u64 {
        self.logical_sample_id
    }

    /// Normalized coordinate inside the shutter interval.
    #[must_use]
    pub const fn normalized_time(self) -> NormalizedShutterTime {
        self.normalized_time
    }

    /// Absolute sample time [s].
    #[must_use]
    pub const fn absolute_time_s(self) -> f64 {
        self.absolute_time_s
    }
}

/// Immutable restart checkpoint containing every bit needed for exact replay.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TemporalAccumulationCheckpoint {
    shutter: ShutterInterval,
    stream_identity: u64,
    pixel_identity: u64,
    color_space: TemporalColorSpace,
    sum: [f64; 3],
    sample_count: u64,
    next_sample_id: u64,
}

impl TemporalAccumulationCheckpoint {
    /// Admitted shutter retained by this accumulation stream.
    #[must_use]
    pub const fn shutter(self) -> ShutterInterval {
        self.shutter
    }

    /// Explicit render stream used to domain-separate shutter samples.
    #[must_use]
    pub const fn stream_identity(self) -> u64 {
        self.stream_identity
    }

    /// Stable pixel/output identity used to domain-separate shutter samples.
    #[must_use]
    pub const fn pixel_identity(self) -> u64 {
        self.pixel_identity
    }

    /// Interpretation of the three accumulated linear channels.
    #[must_use]
    pub const fn color_space(self) -> TemporalColorSpace {
        self.color_space
    }

    /// Number of accepted samples.
    #[must_use]
    pub const fn sample_count(self) -> u64 {
        self.sample_count
    }

    /// Absolute logical sample ID required at the start of the next partition.
    #[must_use]
    pub const fn next_sample_id(self) -> u64 {
        self.next_sample_id
    }

    /// Current channel mean, or `None` before the first accepted sample.
    #[must_use]
    pub fn mean(self) -> Option<[f64; 3]> {
        (self.sample_count != 0).then(|| {
            let denominator = self.sample_count as f64;
            [
                self.sum[0] / denominator,
                self.sum[1] / denominator,
                self.sum[2] / denominator,
            ]
        })
    }
}

/// Mutable accepted state for one temporal color stream.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TemporalAccumulator {
    accepted: TemporalAccumulationCheckpoint,
}

impl TemporalAccumulator {
    /// Start a stream at an explicit absolute logical sample ID.
    #[must_use]
    pub const fn new(
        shutter: ShutterInterval,
        pixel_identity: u64,
        color_space: TemporalColorSpace,
        first_sample_id: u64,
    ) -> Self {
        Self::new_with_stream(shutter, 0, pixel_identity, color_space, first_sample_id)
    }

    /// Start a stream with explicit shutter-stream, pixel, and logical-sample
    /// replay identities.
    #[must_use]
    pub const fn new_with_stream(
        shutter: ShutterInterval,
        stream_identity: u64,
        pixel_identity: u64,
        color_space: TemporalColorSpace,
        first_sample_id: u64,
    ) -> Self {
        Self {
            accepted: TemporalAccumulationCheckpoint {
                shutter,
                stream_identity,
                pixel_identity,
                color_space,
                sum: [0.0; 3],
                sample_count: 0,
                next_sample_id: first_sample_id,
            },
        }
    }

    /// Resume the exact accepted state represented by a typed checkpoint.
    #[must_use]
    pub const fn from_checkpoint(checkpoint: TemporalAccumulationCheckpoint) -> Self {
        Self {
            accepted: checkpoint,
        }
    }

    /// Snapshot every value required for deterministic progressive replay.
    #[must_use]
    pub const fn checkpoint(&self) -> TemporalAccumulationCheckpoint {
        self.accepted
    }

    /// Number of accepted samples.
    #[must_use]
    pub const fn sample_count(&self) -> u64 {
        self.accepted.sample_count
    }

    /// Explicit render stream used to domain-separate shutter samples.
    #[must_use]
    pub const fn stream_identity(&self) -> u64 {
        self.accepted.stream_identity
    }

    /// Absolute logical sample ID required by the next partition.
    #[must_use]
    pub const fn next_sample_id(&self) -> u64 {
        self.accepted.next_sample_id
    }

    /// Current channel mean, or `None` before the first accepted sample.
    #[must_use]
    pub fn mean(&self) -> Option<[f64; 3]> {
        self.accepted.mean()
    }

    /// Evaluate and transactionally append one contiguous absolute-ID range.
    ///
    /// `evaluate` is called exactly once per sample in ascending order.  A
    /// cancellation observed before/after any call, a non-finite sample, an
    /// arithmetic overflow, or a final pre-commit cancellation leaves `self`
    /// bit-for-bit unchanged.
    pub fn accumulate_range(
        &mut self,
        cx: &Cx<'_>,
        sample_ids: Range<u64>,
        mut evaluate: impl FnMut(TemporalSamplePoint) -> [f64; 3],
    ) -> Result<(), TemporalAccumulationError> {
        checkpoint(cx)?;
        if sample_ids.end < sample_ids.start {
            return Err(TemporalAccumulationError::InvalidSampleRange);
        }
        if sample_ids.start != self.accepted.next_sample_id {
            return Err(TemporalAccumulationError::NonContiguousRange {
                expected_start: self.accepted.next_sample_id,
                observed_start: sample_ids.start,
            });
        }
        let additional = sample_ids.end - sample_ids.start;
        let sample_count = self
            .accepted
            .sample_count
            .checked_add(additional)
            .ok_or(TemporalAccumulationError::SampleCountOverflow)?;

        let mut staged = self.accepted;
        for logical_sample_id in sample_ids.clone() {
            checkpoint(cx)?;
            let normalized_time = staged.shutter.sample_for_stream(
                staged.stream_identity,
                staged.pixel_identity,
                logical_sample_id,
            );
            let sample = TemporalSamplePoint {
                logical_sample_id,
                normalized_time,
                absolute_time_s: staged.shutter.time_at(normalized_time),
            };
            let channels = evaluate(sample);
            checkpoint(cx)?;
            for (channel, value) in channels.into_iter().enumerate() {
                if !value.is_finite() {
                    return Err(TemporalAccumulationError::NonFiniteSample {
                        logical_sample_id,
                        channel,
                    });
                }
                let next = staged.sum[channel] + value;
                if !next.is_finite() {
                    return Err(TemporalAccumulationError::AccumulationOverflow {
                        logical_sample_id,
                        channel,
                    });
                }
                staged.sum[channel] = next;
            }
        }
        checkpoint(cx)?;
        staged.sample_count = sample_count;
        staged.next_sample_id = sample_ids.end;
        self.accepted = staged;
        Ok(())
    }
}

/// Typed refusal from one transactional accumulation request.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TemporalAccumulationError {
    /// Cancellation was observed before the transaction committed.
    Cancelled,
    /// `range.end < range.start`.
    InvalidSampleRange,
    /// A partition did not begin at the retained absolute sample checkpoint.
    NonContiguousRange {
        /// Required start ID.
        expected_start: u64,
        /// Supplied start ID.
        observed_start: u64,
    },
    /// The accepted sample counter cannot represent the requested append.
    SampleCountOverflow,
    /// An evaluator returned a NaN or infinity.
    NonFiniteSample {
        /// Absolute logical sample ID whose output was rejected.
        logical_sample_id: u64,
        /// Rejected channel index in `0..3`.
        channel: usize,
    },
    /// A finite sample overflowed the retained finite channel sum.
    AccumulationOverflow {
        /// Absolute logical sample ID whose addition overflowed.
        logical_sample_id: u64,
        /// Overflowed channel index in `0..3`.
        channel: usize,
    },
}

impl fmt::Display for TemporalAccumulationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for TemporalAccumulationError {}

impl From<Cancelled> for TemporalAccumulationError {
    fn from(_: Cancelled) -> Self {
        Self::Cancelled
    }
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), TemporalAccumulationError> {
    cx.checkpoint().map_err(TemporalAccumulationError::from)
}

#[cfg(test)]
mod tests {
    use asupersync::types::Budget;
    use fs_exec::{CancelGate, ExecMode, StreamKey};

    use super::*;
    use crate::motion::{ShotTimeBounds, ShutterConvention, ShutterDistribution};

    fn shutter(duration_s: f64) -> ShutterInterval {
        ShutterInterval::resolve(
            2.0,
            duration_s,
            ShutterConvention::FrontLoaded,
            ShutterDistribution::UniformCounterV1,
            ShotTimeBounds::try_new(0.0, 10.0).expect("shot"),
        )
        .expect("shutter")
    }

    fn with_cx<R>(operation: impl FnOnce(&CancelGate, &Cx<'_>) -> R) -> R {
        let gate = CancelGate::new_clock_free();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 0x5445_4d50,
                    kernel_id: 1,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            operation(&gate, &cx)
        })
    }

    #[test]
    fn g0_constant_rgb_and_xyz_signals_retain_mean_count_and_checkpoint() {
        with_cx(|_, cx| {
            for color_space in [TemporalColorSpace::LinearRgb, TemporalColorSpace::Xyz] {
                let mut accumulator = TemporalAccumulator::new(shutter(1.0), 17, color_space, 40);
                accumulator
                    .accumulate_range(cx, 40..104, |_| [0.25, 2.0, -0.5])
                    .expect("constant signal");
                assert_eq!(accumulator.mean(), Some([0.25, 2.0, -0.5]));
                assert_eq!(accumulator.sample_count(), 64);
                assert_eq!(accumulator.next_sample_id(), 104);
                let checkpoint = accumulator.checkpoint();
                assert_eq!(checkpoint.mean(), accumulator.mean());
                assert_eq!(checkpoint.sample_count(), 64);
                assert_eq!(checkpoint.next_sample_id(), 104);
                assert_eq!(checkpoint.color_space(), color_space);
                assert_eq!(checkpoint.stream_identity(), 0);
                assert_eq!(
                    TemporalAccumulator::from_checkpoint(checkpoint),
                    accumulator
                );
            }
        });
    }

    #[test]
    fn g0_high_rate_linear_time_signal_matches_ordered_reference() {
        with_cx(|_, cx| {
            let shutter = shutter(0.25);
            let mut expected_sum = [0.0; 3];
            for logical_sample_id in 700..4_700 {
                let normalized = shutter.sample(91, logical_sample_id);
                let time = shutter.time_at(normalized);
                let value = [1.0e9 * time - 2.0e9, -7.5e8 * time, 3.0e6 * time + 4.0];
                for channel in 0..3 {
                    expected_sum[channel] += value[channel];
                }
            }
            let expected_mean = expected_sum.map(|sum| sum / 4_000.0);
            let mut accumulator =
                TemporalAccumulator::new(shutter, 91, TemporalColorSpace::Xyz, 700);
            accumulator
                .accumulate_range(cx, 700..4_700, |sample| {
                    let time = sample.absolute_time_s();
                    [1.0e9 * time - 2.0e9, -7.5e8 * time, 3.0e6 * time + 4.0]
                })
                .expect("linear signal");
            assert_eq!(accumulator.mean(), Some(expected_mean));
        });
    }

    #[test]
    fn g5_one_shot_and_progressive_contiguous_partitions_are_bit_identical() {
        with_cx(|_, cx| {
            let shutter = shutter(2.0);
            let evaluate = |sample: TemporalSamplePoint| {
                let id = sample.logical_sample_id() as f64;
                let time = sample.absolute_time_s();
                [time.mul_add(id, 0.125), id.sin(), time.cos()]
            };
            let mut one_shot =
                TemporalAccumulator::new(shutter, 0x55aa, TemporalColorSpace::LinearRgb, 13);
            one_shot
                .accumulate_range(cx, 13..2_013, evaluate)
                .expect("one shot");

            let mut progressive =
                TemporalAccumulator::new(shutter, 0x55aa, TemporalColorSpace::LinearRgb, 13);
            for range in [13..37, 37..911, 911..1_337, 1_337..2_013] {
                progressive
                    .accumulate_range(cx, range, evaluate)
                    .expect("progressive partition");
            }
            assert_eq!(one_shot.checkpoint(), progressive.checkpoint());
            let one_bits = one_shot.mean().expect("mean").map(f64::to_bits);
            let progressive_bits = progressive.mean().expect("mean").map(f64::to_bits);
            assert_eq!(one_bits, progressive_bits);
        });
    }

    #[test]
    fn g0_stream_identity_is_checkpointed_and_reseeds_shutter_coordinates() {
        with_cx(|_, cx| {
            let shutter = shutter(1.0);
            let mut first = TemporalAccumulator::new_with_stream(
                shutter,
                0x1111,
                23,
                TemporalColorSpace::Xyz,
                0,
            );
            let mut replay = first;
            let mut reseeded = TemporalAccumulator::new_with_stream(
                shutter,
                0x2222,
                23,
                TemporalColorSpace::Xyz,
                0,
            );
            let evaluate = |sample: TemporalSamplePoint| {
                let coordinate = sample.normalized_time().value();
                [
                    coordinate,
                    coordinate * coordinate,
                    sample.absolute_time_s(),
                ]
            };
            first.accumulate_range(cx, 0..64, evaluate).unwrap();
            replay.accumulate_range(cx, 0..64, evaluate).unwrap();
            reseeded.accumulate_range(cx, 0..64, evaluate).unwrap();

            assert_eq!(first, replay);
            assert_eq!(first.stream_identity(), 0x1111);
            assert_ne!(first.mean(), reseeded.mean());
        });
    }

    #[test]
    fn g0_pre_and_mid_range_cancellation_leave_state_unchanged() {
        with_cx(|gate, cx| {
            let mut accumulator =
                TemporalAccumulator::new(shutter(1.0), 3, TemporalColorSpace::LinearRgb, 0);
            let original = accumulator;
            gate.request();
            assert_eq!(
                accumulator.accumulate_range(cx, 0..4, |_| [1.0; 3]),
                Err(TemporalAccumulationError::Cancelled)
            );
            assert_eq!(accumulator, original);
        });

        with_cx(|gate, cx| {
            let mut accumulator =
                TemporalAccumulator::new(shutter(1.0), 3, TemporalColorSpace::LinearRgb, 0);
            let original = accumulator;
            let result = accumulator.accumulate_range(cx, 0..8, |sample| {
                if sample.logical_sample_id() == 3 {
                    gate.request();
                }
                [sample.logical_sample_id() as f64; 3]
            });
            assert_eq!(result, Err(TemporalAccumulationError::Cancelled));
            assert_eq!(accumulator, original);
        });
    }

    #[test]
    fn g0_empty_range_is_a_transactional_no_op() {
        with_cx(|_, cx| {
            let mut accumulator =
                TemporalAccumulator::new(shutter(1.0), 5, TemporalColorSpace::Xyz, 99);
            let original = accumulator;
            let mut calls = 0;
            accumulator
                .accumulate_range(cx, 99..99, |_| {
                    calls += 1;
                    [0.0; 3]
                })
                .expect("empty range");
            assert_eq!(calls, 0);
            assert_eq!(accumulator, original);
            assert_eq!(accumulator.mean(), None);
        });
    }

    #[test]
    fn g0_overflow_and_non_finite_output_roll_back_the_whole_range() {
        with_cx(|_, cx| {
            let mut overflow =
                TemporalAccumulator::new(shutter(1.0), 8, TemporalColorSpace::LinearRgb, 0);
            overflow
                .accumulate_range(cx, 0..1, |_| [f64::MAX, 0.0, 0.0])
                .expect("first finite maximum");
            let accepted = overflow;
            assert_eq!(
                overflow.accumulate_range(cx, 1..3, |_| [f64::MAX, 1.0, 1.0]),
                Err(TemporalAccumulationError::AccumulationOverflow {
                    logical_sample_id: 1,
                    channel: 0
                })
            );
            assert_eq!(overflow, accepted);

            let mut malformed =
                TemporalAccumulator::new(shutter(1.0), 8, TemporalColorSpace::Xyz, 10);
            let original = malformed;
            assert_eq!(
                malformed.accumulate_range(cx, 10..14, |sample| {
                    if sample.logical_sample_id() == 12 {
                        [1.0, f64::NAN, 3.0]
                    } else {
                        [1.0, 2.0, 3.0]
                    }
                }),
                Err(TemporalAccumulationError::NonFiniteSample {
                    logical_sample_id: 12,
                    channel: 1
                })
            );
            assert_eq!(malformed, original);
        });
    }

    #[test]
    #[allow(clippy::reversed_empty_ranges)] // deliberately malformed public input
    fn g0_range_and_count_overflow_refuse_without_mutation() {
        with_cx(|_, cx| {
            let mut accumulator =
                TemporalAccumulator::new(shutter(1.0), 1, TemporalColorSpace::Xyz, 7);
            let original = accumulator;
            assert_eq!(
                accumulator.accumulate_range(cx, 8..7, |_| [0.0; 3]),
                Err(TemporalAccumulationError::InvalidSampleRange)
            );
            assert_eq!(accumulator, original);
            assert_eq!(
                accumulator.accumulate_range(cx, 8..9, |_| [0.0; 3]),
                Err(TemporalAccumulationError::NonContiguousRange {
                    expected_start: 7,
                    observed_start: 8
                })
            );
            assert_eq!(accumulator, original);

            accumulator.accepted.sample_count = u64::MAX;
            let saturated = accumulator;
            assert_eq!(
                accumulator.accumulate_range(cx, 7..8, |_| [0.0; 3]),
                Err(TemporalAccumulationError::SampleCountOverflow)
            );
            assert_eq!(accumulator, saturated);
        });
    }

    #[test]
    fn g0_zero_width_shutter_evaluates_every_id_at_one_exact_time() {
        with_cx(|_, cx| {
            let shutter = shutter(0.0);
            let exact_time = shutter.open_s();
            let mut accumulator =
                TemporalAccumulator::new(shutter, 44, TemporalColorSpace::LinearRgb, 100);
            accumulator
                .accumulate_range(cx, 100..132, |sample| {
                    assert_eq!(sample.absolute_time_s().to_bits(), exact_time.to_bits());
                    [sample.absolute_time_s(); 3]
                })
                .expect("zero-width shutter");
            assert_eq!(accumulator.mean(), Some([exact_time; 3]));
        });
    }
}
