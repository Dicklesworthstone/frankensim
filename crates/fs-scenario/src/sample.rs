//! Shared time sampling for scalar histories and immutable payload carriers.

use crate::ScenarioError;
use crate::payload::{OutsideDomainPolicy, SampleSource, TableInterpolation};
use crate::signal::TimeSignal;
use fs_qty::QtyAny;

/// Admission limits for a uniformly sampled time history.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SamplingBudget {
    /// Maximum number of returned samples.
    pub max_samples: usize,
    /// Maximum storage for the returned `QtyAny` buffer, including dimensions.
    /// The borrowed source and bounded validation diagnostics are not counted.
    pub max_workspace_bytes: usize,
}

/// Default uniform-history storage admission: at most one million samples.
pub const DEFAULT_SAMPLING_BUDGET: SamplingBudget = SamplingBudget {
    max_samples: 1_000_000,
    max_workspace_bytes: 64 * 1024 * 1024,
};

impl Default for SamplingBudget {
    fn default() -> Self {
        DEFAULT_SAMPLING_BUDGET
    }
}

/// A refused uniform history never publishes a partial sample vector.
#[derive(Debug, Clone, PartialEq)]
pub enum SamplingError {
    /// Nonfinite, nonpositive, overflowing or non-advancing sample clock.
    InvalidClock(&'static str),
    /// A declared source frequency reaches or exceeds half the sample rate.
    Nyquist {
        /// Largest instantaneous source frequency on the sampled interval.
        maximum_frequency_hz: f64,
        /// Half the requested sample rate.
        nyquist_hz: f64,
    },
    /// The requested sample count or buffer size exceeds its explicit limit.
    Resource {
        /// Refused resource: `samples` or `workspace_bytes`.
        resource: &'static str,
        /// Required amount, computed without overflowing machine-size counts.
        requested: u128,
        /// Caller-supplied limit for the resource.
        limit: u128,
    },
    /// The admitted sample buffer could not be allocated.
    Allocation,
    /// The original signal owner refused validation or evaluation.
    Signal(ScenarioError),
}

impl std::fmt::Display for SamplingError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "scenario uniform sampling refused: {self:?}")
    }
}

impl std::error::Error for SamplingError {}

impl TimeSignal {
    /// Sample this time history on an explicitly admitted uniform clock.
    ///
    /// Samples use `start_s + index / sample_rate_hz`, with separate division
    /// and addition and no accumulated time step. The source's phase and
    /// physical frequencies remain independent of this observation clock.
    /// A chirp's maximum instantaneous frequency over the first-to-last sample
    /// interval must be strictly below Nyquist. This is source-frequency
    /// admission, not a claim that a finite chirp record is bandlimited or that
    /// unspecified table/Chebyshev histories have been antialias filtered.
    ///
    /// Zero samples validates the clock and budget, then returns an empty
    /// vector without inspecting or evaluating the signal. A nonempty request
    /// validates the signal once before sampling through its original owner.
    ///
    /// # Errors
    /// Invalid clocks, exhausted budgets, allocation failure, source frequency
    /// at or above Nyquist, and ordinary signal validation/evaluation errors.
    pub fn sample_uniform(
        &self,
        start_s: f64,
        sample_rate_hz: f64,
        sample_count: usize,
        budget: SamplingBudget,
    ) -> Result<Vec<QtyAny>, SamplingError> {
        match self.sample_uniform_with_checkpoint(
            start_s,
            sample_rate_hz,
            sample_count,
            budget,
            &mut |_| Ok::<(), core::convert::Infallible>(()),
        ) {
            Ok(result) => result,
            Err(never) => match never {},
        }
    }

    /// Cancellable uniform sampling with the same clock and admission rules as
    /// [`Self::sample_uniform`]. The outer result retains the caller's exact
    /// checkpoint error. Validation, every sample, spectral recurrences and
    /// final publication are cancellable; no refusal returns a partial vector.
    pub fn sample_uniform_with_checkpoint<E>(
        &self,
        start_s: f64,
        sample_rate_hz: f64,
        sample_count: usize,
        budget: SamplingBudget,
        checkpoint: &mut impl FnMut(&'static str) -> Result<(), E>,
    ) -> Result<Result<Vec<QtyAny>, SamplingError>, E> {
        checkpoint("uniform sampling admission")?;
        if !start_s.is_finite() || !sample_rate_hz.is_finite() || sample_rate_hz <= 0.0 {
            return Ok(Err(SamplingError::InvalidClock(
                "finite start time and positive finite sample rate required",
            )));
        }
        if sample_count > budget.max_samples {
            return Ok(Err(SamplingError::Resource {
                resource: "samples",
                requested: sample_count as u128,
                limit: budget.max_samples as u128,
            }));
        }
        let bytes = (sample_count as u128) * (std::mem::size_of::<QtyAny>() as u128);
        if bytes > budget.max_workspace_bytes as u128 {
            return Ok(Err(SamplingError::Resource {
                resource: "workspace_bytes",
                requested: bytes,
                limit: budget.max_workspace_bytes as u128,
            }));
        }
        if sample_count == 0 {
            checkpoint("uniform sampling publication")?;
            return Ok(Ok(Vec::new()));
        }
        if (sample_count - 1) as u128 > (1_u128 << 53) {
            return Ok(Err(SamplingError::InvalidClock(
                "sample index is not exactly representable on the clock",
            )));
        }
        let end_s = start_s + (sample_count - 1) as f64 / sample_rate_hz;
        if !end_s.is_finite()
            || (sample_count > 1
                && (start_s + 1.0 / sample_rate_hz <= start_s
                    || start_s + (sample_count - 2) as f64 / sample_rate_hz >= end_s))
        {
            return Ok(Err(SamplingError::InvalidClock(
                "sample times must remain finite and strictly advance",
            )));
        }
        if let Self::Chirp(chirp) = self {
            let maximum_frequency_hz = match chirp.max_frequency_hz(start_s, end_s) {
                Ok(frequency) => frequency,
                Err(error) => return Ok(Err(SamplingError::Signal(error))),
            };
            // The ratio also handles a positive subnormal rate whose half
            // underflows to zero; a DC source still lies below that Nyquist.
            if maximum_frequency_hz / sample_rate_hz >= 0.5 {
                return Ok(Err(SamplingError::Nyquist {
                    maximum_frequency_hz,
                    nyquist_hz: sample_rate_hz * 0.5,
                }));
            }
        }
        let mut violations = Vec::new();
        self.check_with_checkpoint("uniform time sampling", &mut violations, checkpoint)?;
        if let Some(violation) = violations.into_iter().next() {
            return Ok(Err(SamplingError::Signal(ScenarioError::Evaluate {
                what: violation.what,
            })));
        }
        checkpoint("uniform sampling allocation")?;
        let mut samples = Vec::new();
        if samples.try_reserve_exact(sample_count).is_err() {
            return Ok(Err(SamplingError::Allocation));
        }
        let mut previous_time = None;
        for index in 0..sample_count {
            checkpoint("uniform sampling sample")?;
            let time_s = start_s + index as f64 / sample_rate_hz;
            if previous_time.is_some_and(|previous| time_s <= previous) {
                return Ok(Err(SamplingError::InvalidClock(
                    "sample times must remain finite and strictly advance",
                )));
            }
            let sample = match self.eval_prevalidated_with_checkpoint(time_s, &mut || {
                checkpoint("uniform sampling Chebyshev evaluation")
            })? {
                Ok(sample) => sample,
                Err(error) => return Ok(Err(SamplingError::Signal(error))),
            };
            samples.push(sample);
            previous_time = Some(time_s);
        }
        checkpoint("uniform sampling publication")?;
        Ok(Ok(samples))
    }
}

/// No sample is returned after a policy, arithmetic, frame or caller refusal.
#[derive(Debug, Clone, PartialEq)]
pub enum PayloadEvaluationError {
    Unsupported(&'static str),
    InvalidTime,
    OutsideDomain {
        time_s: f64,
        first_s: f64,
        last_s: f64,
    },
    Arithmetic(&'static str),
    Frame(ScenarioError),
    Cancelled,
}
impl std::fmt::Display for PayloadEvaluationError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "scenario payload evaluation refused: {self:?}")
    }
}
impl std::error::Error for PayloadEvaluationError {}

pub(crate) fn value<T, const N: usize>(
    source: &SampleSource<T>,
    time_s: f64,
    components: impl Fn(&T) -> [f64; N],
) -> Result<[f64; N], PayloadEvaluationError> {
    match source {
        SampleSource::Fixed(value) => Ok(components(value)),
        SampleSource::Distribution { .. } => Err(PayloadEvaluationError::Unsupported(
            "a distribution is not an explicit time history",
        )),
        SampleSource::Table {
            times,
            values,
            interpolation,
            outside,
        } => table_value(times, values, *interpolation, *outside, time_s, components),
    }
}

fn table_value<T, const N: usize>(
    times: &[QtyAny],
    values: &[T],
    interpolation: TableInterpolation,
    outside: OutsideDomainPolicy,
    time_s: f64,
    components: impl Fn(&T) -> [f64; N],
) -> Result<[f64; N], PayloadEvaluationError> {
    use PayloadEvaluationError as E;
    // Nonempty, increasing, finite time and equal shapes are constructor invariants.
    let first = times[0].value;
    let last = times[times.len() - 1].value;
    let mut time = time_s;
    if time < first || time > last {
        time = match outside {
            OutsideDomainPolicy::Refuse => {
                return Err(E::OutsideDomain {
                    time_s,
                    first_s: first,
                    last_s: last,
                });
            }
            OutsideDomainPolicy::Clamp => time.clamp(first, last),
            OutsideDomainPolicy::Periodic => {
                let period = last - first;
                let delta = time - first;
                if !period.is_finite() || period <= 0.0 || !delta.is_finite() {
                    return Err(E::Arithmetic("unrepresentable periodic time reduction"));
                }
                let phase = delta.rem_euclid(period);
                let wrapped = first + phase;
                if !wrapped.is_finite() {
                    return Err(E::Arithmetic("nonfinite wrapped time"));
                }
                wrapped.clamp(first, last)
            }
        };
    }
    let left = times.partition_point(|t| t.value <= time).saturating_sub(1);
    if left + 1 == times.len()
        || time == times[left].value
        || interpolation == TableInterpolation::StepLeft
    {
        return Ok(components(&values[left]));
    }
    let width = times[left + 1].value - times[left].value;
    if !width.is_finite() || width <= 0.0 {
        return Err(E::Arithmetic("unrepresentable interpolation span"));
    }
    let weight = (time - times[left].value) / width;
    if !weight.is_finite() || !(0.0..=1.0).contains(&weight) {
        return Err(E::Arithmetic("unrepresentable interpolation weight"));
    }
    let a = components(&values[left]);
    let b = components(&values[left + 1]);
    // Preserve the original vector evaluator's convex arithmetic and exact knots.
    let out = std::array::from_fn(|i| (1.0 - weight) * a[i] + weight * b[i]);
    if out.iter().all(|v| v.is_finite()) {
        Ok(out)
    } else {
        Err(E::Arithmetic("nonfinite interpolated component"))
    }
}

#[cfg(test)]
mod uniform_tests {
    use super::{DEFAULT_SAMPLING_BUDGET, SamplingBudget, SamplingError};
    use crate::signal::{Interp, LinearChirp, TimeSignal};
    use fs_qty::{Dims, QtyAny};

    fn chirp(from_hz: f64, to_hz: f64) -> TimeSignal {
        TimeSignal::Chirp(
            LinearChirp::new(
                0.0,
                1.0,
                from_hz,
                to_hz,
                0.0,
                QtyAny::new(2.0, Dims([1, 1, -2, 0, 0, 0])),
                QtyAny::new(3.0, Dims([1, 1, -2, 0, 0, 0])),
            )
            .unwrap(),
        )
    }

    #[test]
    fn g0_frequency_admission_covers_the_actual_sweep_interval_and_tails() {
        let rising = chirp(2.0, 10.0);
        // A low start frequency cannot conceal a later Nyquist crossing.
        assert_eq!(
            rising.sample_uniform(0.0, 16.0, 13, DEFAULT_SAMPLING_BUDGET),
            Err(SamplingError::Nyquist {
                maximum_frequency_hz: 8.0,
                nyquist_hz: 8.0,
            })
        );
        // The source's later, unsampled high frequencies do not reject a prefix.
        assert_eq!(
            rising
                .sample_uniform(0.0, 16.0, 12, DEFAULT_SAMPLING_BUDGET)
                .unwrap()
                .len(),
            12
        );
        assert!(matches!(
            rising.sample_uniform(2.0, 16.0, 1, DEFAULT_SAMPLING_BUDGET),
            Err(SamplingError::Nyquist { .. })
        ));
        let falling = chirp(10.0, 2.0);
        assert!(
            falling
                .sample_uniform(0.5, 16.0, 17, DEFAULT_SAMPLING_BUDGET)
                .is_ok()
        );
        assert!(matches!(
            falling.sample_uniform(-1.0, 16.0, 1, DEFAULT_SAMPLING_BUDGET),
            Err(SamplingError::Nyquist { .. })
        ));
        assert!(
            chirp(f64::from_bits(8.0_f64.to_bits() - 1), 0.0)
                .sample_uniform(0.0, 16.0, 1, DEFAULT_SAMPLING_BUDGET)
                .is_ok()
        );
    }

    #[test]
    fn g3_clock_refinement_preserves_physical_phase_values_and_dimensions() {
        for signal in [chirp(2.0, 2.0), chirp(2.0, 6.0)] {
            let coarse = signal
                .sample_uniform(-0.25, 16.0, 25, DEFAULT_SAMPLING_BUDGET)
                .unwrap();
            let fine = signal
                .sample_uniform(-0.25, 32.0, 49, DEFAULT_SAMPLING_BUDGET)
                .unwrap();
            assert_eq!(coarse[4].value, 5.0, "phase zero is the positive peak");
            for (index, value) in coarse.iter().enumerate() {
                assert_eq!(*value, fine[2 * index]);
                assert_eq!(*value, signal.eval(-0.25 + index as f64 / 16.0).unwrap());
                assert_eq!(value.dims, Dims([1, 1, -2, 0, 0, 0]));
            }
        }
        let tone = chirp(2.0, 2.0)
            .sample_uniform(0.0, 16.0, 9, DEFAULT_SAMPLING_BUDGET)
            .unwrap();
        assert_eq!(tone[0].value, 5.0);
        assert!((tone[2].value - 2.0).abs() < 1e-14);
        assert_eq!(tone[4].value, -1.0);
        assert_eq!(tone[8].value, 5.0);
    }

    #[test]
    fn g0_invalid_and_nonadvancing_clocks_are_refused_before_allocation() {
        let signal = TimeSignal::Constant(QtyAny::dimensionless(1.0));
        for (start, rate, count) in [
            (f64::NAN, 16.0, 1),
            (f64::INFINITY, 16.0, 1),
            (0.0, 0.0, 1),
            (0.0, -16.0, 1),
            (0.0, f64::NAN, 1),
            (0.0, f64::INFINITY, 1),
            (9_007_199_254_740_992.0, 1.0, 2),
            (f64::MAX, f64::MIN_POSITIVE, 3),
        ] {
            assert!(matches!(
                signal.sample_uniform(start, rate, count, DEFAULT_SAMPLING_BUDGET),
                Err(SamplingError::InvalidClock(_))
            ));
        }
        // Positive subnormal clocks must not misclassify DC as Nyquist.
        assert_eq!(
            chirp(0.0, 0.0)
                .sample_uniform(0.0, f64::from_bits(1), 1, DEFAULT_SAMPLING_BUDGET)
                .unwrap()[0]
                .value,
            5.0
        );
    }

    #[test]
    fn g0_sample_and_full_quantity_storage_budgets_are_independent() {
        let signal = chirp(1.0, 1.0);
        let budget = SamplingBudget {
            max_samples: 1,
            ..DEFAULT_SAMPLING_BUDGET
        };
        assert!(matches!(
            signal.sample_uniform(0.0, 16.0, 2, budget),
            Err(SamplingError::Resource {
                resource: "samples",
                requested: 2,
                limit: 1
            })
        ));
        let bytes = 2 * std::mem::size_of::<QtyAny>();
        let budget = SamplingBudget {
            max_samples: 2,
            max_workspace_bytes: bytes - 1,
        };
        assert!(matches!(
            signal.sample_uniform(0.0, 16.0, 2, budget),
            Err(SamplingError::Resource {
                resource: "workspace_bytes",
                ..
            })
        ));
        assert_eq!(
            signal
                .sample_uniform(
                    0.0,
                    16.0,
                    2,
                    SamplingBudget {
                        max_workspace_bytes: bytes,
                        ..budget
                    }
                )
                .unwrap()
                .len(),
            2
        );
        assert!(matches!(
            signal.sample_uniform(
                0.0,
                16.0,
                usize::MAX,
                SamplingBudget {
                    max_samples: usize::MAX,
                    max_workspace_bytes: usize::MAX
                }
            ),
            Err(SamplingError::Resource {
                resource: "workspace_bytes",
                ..
            })
        ));
    }

    #[test]
    fn g0_table_history_is_validated_once_and_uses_original_knot_semantics() {
        let signal = TimeSignal::Table {
            times: vec![0.0, 0.5, 1.0],
            values: vec![1.0, 2.0, 3.0],
            dims: Dims::NONE,
            interp: Interp::Hold,
        };
        let mut validation_visits = 0;
        let samples = signal
            .sample_uniform_with_checkpoint(-0.25, 4.0, 7, DEFAULT_SAMPLING_BUDGET, &mut |phase| {
                validation_visits += usize::from(matches!(
                    phase,
                    "signal table times" | "signal table values"
                ));
                Ok::<(), ()>(())
            })
            .unwrap()
            .unwrap();
        assert_eq!(validation_visits, 6);
        assert_eq!(
            samples.into_iter().map(|q| q.value).collect::<Vec<_>>(),
            vec![1.0, 1.0, 1.0, 2.0, 2.0, 3.0, 3.0]
        );
        let mut malformed = signal;
        if let TimeSignal::Table { values, .. } = &mut malformed {
            values[2] = f64::NAN;
        }
        assert!(matches!(
            malformed.sample_uniform(-1.0, 4.0, 1, DEFAULT_SAMPLING_BUDGET),
            Err(SamplingError::Signal(_))
        ));
    }

    #[test]
    fn g4_cancellation_never_returns_a_partial_sample_set() {
        for stop_phase in [
            "uniform sampling admission",
            "uniform sampling sample",
            "uniform sampling publication",
        ] {
            let result = chirp(2.0, 4.0).sample_uniform_with_checkpoint(
                0.0,
                16.0,
                8,
                DEFAULT_SAMPLING_BUDGET,
                &mut |phase| {
                    if phase == stop_phase {
                        Err("caller cancelled")
                    } else {
                        Ok(())
                    }
                },
            );
            assert_eq!(result, Err("caller cancelled"));
        }
        let table = TimeSignal::Table {
            times: vec![0.0, 1.0],
            values: vec![1.0, 2.0],
            dims: Dims::NONE,
            interp: Interp::Linear,
        };
        let result = table.sample_uniform_with_checkpoint(
            0.0,
            16.0,
            8,
            DEFAULT_SAMPLING_BUDGET,
            &mut |phase| {
                if phase == "signal table values" {
                    Err(23)
                } else {
                    Ok(())
                }
            },
        );
        assert_eq!(result, Err(23));
    }

    #[test]
    fn g0_empty_requests_validate_only_clock_and_budget() {
        let malformed = TimeSignal::Constant(QtyAny::dimensionless(f64::NAN));
        let no_storage = SamplingBudget {
            max_samples: 0,
            max_workspace_bytes: 0,
        };
        assert!(
            malformed
                .sample_uniform(0.0, 16.0, 0, no_storage)
                .unwrap()
                .is_empty()
        );
        assert!(
            chirp(200.0, 400.0)
                .sample_uniform(0.0, 16.0, 0, no_storage)
                .unwrap()
                .is_empty()
        );
        assert!(matches!(
            malformed.sample_uniform(0.0, 0.0, 0, no_storage),
            Err(SamplingError::InvalidClock(_))
        ));
    }
}
