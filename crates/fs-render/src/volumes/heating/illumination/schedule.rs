//! Finite, piecewise-linear incident power on the native physical clock.
//!
//! Powers are already wavelength/quadrature-weighted watts. Segments may meet
//! with a jump or leave a dark gap; power is exactly zero outside their support.
//! Interval means integrate local overlaps, never subtract cumulative fluences,
//! so a short late pulse is not lost behind a much larger earlier exposure.

use super::super::feedback::OpticalFeedbackError;
use super::super::{CheckedSum, nonnegative};
use fs_blake3::{ContentHash, DomainHasher};
use fs_thermal::enthalpy::DepositionError;
use std::sync::Arc;

/// One linear power segment, with times in seconds and independent lane watts.
/// Endpoint jumps between adjoining segments have zero integrated duration.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpectralPowerSegment<const N: usize> {
    pub start_s: f64,
    pub end_s: f64,
    pub start_power_w: [f64; N],
    pub end_power_w: [f64; N],
}

#[derive(Debug)]
struct Model<const N: usize> {
    identity: ContentHash,
    segments: Vec<SpectralPowerSegment<N>>,
    peak_power_w: [f64; N],
}

/// Immutable, bounded source data. Clones share the original segment storage.
/// This is a prescribed waveform, not an inferred laser or material model.
#[derive(Clone, Debug)]
pub struct SpectralPowerSchedule<const N: usize>(Arc<Model<N>>);

impl<const N: usize> SpectralPowerSchedule<N> {
    /// Admit ordered, nonoverlapping segments. An empty schedule is all-dark.
    pub fn try_new(
        segments: &[SpectralPowerSegment<N>],
        max_segments: usize,
    ) -> Result<Self, OpticalFeedbackError> {
        if segments.len() > max_segments {
            return Err(OpticalFeedbackError::Budget);
        }
        if N == 0
            || segments.iter().any(|s| {
                !nonnegative(s.start_s)
                    || !nonnegative(s.end_s)
                    || s.end_s <= s.start_s
                    || !s
                        .start_power_w
                        .iter()
                        .chain(&s.end_power_w)
                        .copied()
                        .all(nonnegative)
            })
            || segments.windows(2).any(|s| s[0].end_s > s[1].start_s)
        {
            return Err(invalid(
                "ordered nonoverlapping positive-duration segments and finite nonnegative seconds/watts required",
            ));
        }
        let mut retained = Vec::new();
        retained
            .try_reserve_exact(segments.len())
            .map_err(|_| OpticalFeedbackError::Budget)?;
        let mut peak_power_w = [0.0_f64; N];
        let mut hash = DomainHasher::new("org.frankensim.render.spectral-power-schedule.v1");
        hash.update(&(N as u64).to_le_bytes());
        hash.update(&(segments.len() as u64).to_le_bytes());
        for segment in segments {
            for value in [segment.start_s, segment.end_s]
                .into_iter()
                .chain(segment.start_power_w)
                .chain(segment.end_power_w)
            {
                hash.update(&value.to_bits().to_le_bytes());
            }
            for (lane, peak) in peak_power_w.iter_mut().enumerate() {
                *peak = peak
                    .max(segment.start_power_w[lane])
                    .max(segment.end_power_w[lane]);
            }
            retained.push(*segment);
        }
        Ok(Self(Arc::new(Model {
            identity: hash.finalize(),
            segments: retained,
            peak_power_w,
        })))
    }

    pub fn identity(&self) -> ContentHash {
        self.0.identity
    }

    /// Per-lane envelope for transport admission; peaks need not be simultaneous.
    pub fn peak_power_w(&self) -> [f64; N] {
        self.0.peak_power_w
    }

    pub fn segments(&self) -> &[SpectralPowerSegment<N>] {
        &self.0.segments
    }

    pub fn segment_count(&self) -> usize {
        self.0.segments.len()
    }

    /// First support/linear-segment boundary strictly after `after_s`.
    /// Lookup is logarithmic; adjoining jumps share one transition instant.
    pub fn next_transition_time_s(
        &self,
        after_s: f64,
    ) -> Result<Option<f64>, OpticalFeedbackError> {
        if !nonnegative(after_s) {
            return Err(invalid("finite nonnegative source query time required"));
        }
        let index = self.0.segments.partition_point(|s| s.end_s <= after_s);
        Ok(self.0.segments.get(index).map(|s| {
            if s.start_s > after_s {
                s.start_s
            } else {
                s.end_s
            }
        }))
    }

    /// Exact linear-profile mean up to floating-point arithmetic on the
    /// represented `[start_s, end_s]`. Multiply by the original native duration
    /// for energy: clock rounding can make `end_s - start_s` differ from it.
    /// Work is a binary lookup plus the overlapping segments, with cancellation
    /// at each segment. There is no accumulated-fluence state to drift on resume.
    pub fn mean_power_w(
        &self,
        start_s: f64,
        duration_s: f64,
        end_s: f64,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<[f64; N], OpticalFeedbackError> {
        checkpoint()?;
        if !nonnegative(start_s)
            || !nonnegative(end_s)
            || !duration_s.is_finite()
            || duration_s <= 0.0
            || end_s <= start_s
            || (start_s + duration_s != end_s && duration_s != end_s - start_s)
        {
            return Err(invalid("matching finite native source interval required"));
        }
        let width = end_s - start_s;
        let first = self.0.segments.partition_point(|s| s.end_s <= start_s);
        let mut sums = [CheckedSum::default(); N];
        for segment in &self.0.segments[first..] {
            if segment.start_s >= end_s {
                break;
            }
            checkpoint()?;
            let left = start_s.max(segment.start_s);
            let right = end_s.min(segment.end_s);
            let span = segment.end_s - segment.start_s;
            let fraction = (right - left) / width;
            let u = (left - segment.start_s) / span;
            let v = (right - segment.start_s) / span;
            for (lane, sum) in sums.iter_mut().enumerate() {
                let a = interpolate(segment.start_power_w[lane], segment.end_power_w[lane], u);
                let b = interpolate(segment.start_power_w[lane], segment.end_power_w[lane], v);
                // A convex mean avoids overflow of a+b and cancellation of a-b.
                let mean = a.min(b) + 0.5 * (a.max(b) - a.min(b));
                let contribution = mean * fraction;
                if !nonnegative(contribution) || (mean > 0.0 && contribution == 0.0) {
                    return Err(invalid("unrepresentable interval source power"));
                }
                sum.add(contribution)?;
            }
        }
        let mut result = [0.0; N];
        for lane in 0..N {
            let value = sums[lane].value()?;
            let peak = self.0.peak_power_w[lane];
            // Only roundoff of the compensated convex combination may exceed
            // its envelope. Do not hide a physically larger source by clipping.
            if value > peak && value - peak > 8.0 * f64::EPSILON * peak {
                return Err(invalid("interval mean exceeds admitted source envelope"));
            }
            result[lane] = value.min(peak);
        }
        checkpoint()?;
        Ok(result)
    }
}

fn interpolate(a: f64, b: f64, fraction: f64) -> f64 {
    if fraction == 0.0 {
        a
    } else if fraction == 1.0 {
        b
    } else if a <= b {
        a + (b - a) * fraction
    } else {
        b + (a - b) * (1.0 - fraction)
    }
}

fn invalid(message: &'static str) -> OpticalFeedbackError {
    OpticalFeedbackError::Invalid(message)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn independent_spectral_fluences_and_jump_boundaries() {
        let profile = SpectralPowerSchedule::try_new(
            &[
                SpectralPowerSegment {
                    start_s: 1.0,
                    end_s: 2.0,
                    start_power_w: [4.0, 0.0],
                    end_power_w: [4.0, 2.0],
                },
                SpectralPowerSegment {
                    start_s: 2.0,
                    end_s: 3.0,
                    start_power_w: [0.0, 8.0],
                    end_power_w: [6.0, 0.0],
                },
                SpectralPowerSegment {
                    start_s: 4.0,
                    end_s: 5.0,
                    start_power_w: [2.0, 0.0],
                    end_power_w: [0.0, 0.0],
                },
            ],
            3,
        )
        .unwrap();
        assert_eq!(profile.peak_power_w(), [6.0, 8.0]);
        assert_eq!(
            profile.mean_power_w(0.0, 5.0, 5.0, || Ok(())).unwrap(),
            [1.6, 1.0]
        );
        assert_eq!(
            profile.mean_power_w(1.5, 1.0, 2.5, || Ok(())).unwrap(),
            [2.75, 3.75]
        );
        assert_eq!(
            profile.mean_power_w(3.0, 1.0, 4.0, || Ok(())).unwrap(),
            [0.0; 2]
        );
        for (time, expected) in [
            (0.0, Some(1.0)),
            (1.0, Some(2.0)),
            (2.0, Some(3.0)),
            (3.0, Some(4.0)),
            (5.0, None),
        ] {
            assert_eq!(profile.next_transition_time_s(time).unwrap(), expected);
        }
    }

    #[test]
    fn tiny_late_pulse_is_independent_of_huge_earlier_fluence() {
        let late = 2.0_f64.powi(40);
        let pulse = SpectralPowerSegment {
            start_s: late,
            end_s: late + 0.5,
            start_power_w: [2.0e-12],
            end_power_w: [6.0e-12],
        };
        let profile = SpectralPowerSchedule::try_new(
            &[
                SpectralPowerSegment {
                    start_s: 0.0,
                    end_s: 1.0,
                    start_power_w: [1.0e200],
                    end_power_w: [1.0e200],
                },
                pulse,
            ],
            2,
        )
        .unwrap();
        let isolated = SpectralPowerSchedule::try_new(&[pulse], 1).unwrap();
        assert_eq!(
            profile
                .mean_power_w(late, 0.5, late + 0.5, || Ok(()))
                .unwrap(),
            isolated
                .mean_power_w(late, 0.5, late + 0.5, || Ok(()))
                .unwrap()
        );
        let mean = profile
            .mean_power_w(late, 0.5, late + 0.5, || Ok(()))
            .unwrap()[0];
        assert!((mean - 4.0e-12).abs() <= 2.0 * f64::EPSILON * 4.0e-12);
        // The source means live on the represented clock span, even when the
        // native duration is not exactly recoverable by subtracting endpoints.
        let dt = 0.1;
        assert_ne!((late + dt) - late, dt);
        assert!(profile.mean_power_w(late, dt, late + dt, || Ok(())).is_ok());
    }

    #[test]
    fn bounded_admission_cancellation_and_dark_schedule() {
        let segment = SpectralPowerSegment {
            start_s: 0.0,
            end_s: 1.0,
            start_power_w: [1.0],
            end_power_w: [2.0],
        };
        assert!(matches!(
            SpectralPowerSchedule::try_new(&[segment], 0),
            Err(OpticalFeedbackError::Budget)
        ));
        assert!(SpectralPowerSchedule::try_new(&[segment, segment], 2).is_err());
        assert!(
            SpectralPowerSchedule::try_new(
                &[SpectralPowerSegment {
                    end_power_w: [f64::NAN],
                    ..segment
                }],
                1
            )
            .is_err()
        );
        let source = SpectralPowerSchedule::try_new(&[segment], 1).unwrap();
        let mut calls = 0;
        assert!(matches!(
            source.mean_power_w(0.0, 1.0, 1.0, || {
                calls += 1;
                if calls == 2 {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            }),
            Err(OpticalFeedbackError::Thermal(DepositionError::Cancelled))
        ));
        assert!(source.mean_power_w(0.0, 0.5, 1.0, || Ok(())).is_err());
        assert!(source.next_transition_time_s(f64::INFINITY).is_err());
        let dark = SpectralPowerSchedule::<2>::try_new(&[], 0).unwrap();
        assert_eq!(
            dark.mean_power_w(0.0, 1.0, 1.0, || Ok(())).unwrap(),
            [0.0; 2]
        );
        assert_eq!(dark.next_transition_time_s(0.0).unwrap(), None);
        assert_eq!(source.identity(), source.clone().identity());
        assert_ne!(
            source.identity(),
            SpectralPowerSchedule::try_new(
                &[SpectralPowerSegment {
                    end_power_w: [3.0],
                    ..segment
                }],
                1
            )
            .unwrap()
            .identity()
        );
    }
}
