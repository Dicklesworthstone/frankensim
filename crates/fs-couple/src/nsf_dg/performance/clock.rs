//! Conservative overlap of two uniform physical time partitions.
//! Endpoints are integer-index products in output-interval units, never a
//! repeatedly accumulated floating-point clock. Values on the source partition
//! are its original interval means, held constant only for this projection.

pub(crate) const EXACT_INTEGER_LIMIT: u64 = 1 << 53;

#[derive(Clone, Copy, Debug)]
pub(crate) struct IntervalClock {
    width: f64,
    source_intervals: usize,
    output_samples: u64,
}

impl IntervalClock {
    pub(crate) fn new(dt_s: f64, output_rate: u32, output_samples: u64)
        -> Result<Self, &'static str>
    {
        let width = dt_s * f64::from(output_rate);
        // Every integer endpoint must be representable. Real retained-source
        // memory/work budgets normally impose much smaller limits.
        if !dt_s.is_finite() || dt_s <= 0.0 || output_rate == 0
            || !width.is_finite() || width <= 0.0 || output_samples == 0
            || output_samples > EXACT_INTEGER_LIMIT
        {
            return Err("resolved observer needs positive finite clocks and exactly representable sample indices");
        }
        let end = output_samples as f64;
        let count = (end / width).ceil();
        if !count.is_finite() || count < 1.0 || count > EXACT_INTEGER_LIMIT as f64 {
            return Err("resolved observer physical interval count is not representable");
        }
        let mut count = count as u64;
        // Division can round onto an integer on either side. Admit coverage
        // using the SAME endpoint products later used by overlap quadrature.
        if count > 1 && (count - 1) as f64 * width >= end { count -= 1; }
        if count as f64 * width < end {
            count = count.checked_add(1).ok_or("resolved physical interval overflow")?;
        }
        if count > EXACT_INTEGER_LIMIT || count as f64 * width < end
            || (count > 1 && (count - 1) as f64 * width >= end)
        {
            return Err("resolved observer cannot distinguish adjacent physical interval endpoints");
        }
        let source_intervals = usize::try_from(count)
            .map_err(|_| "resolved physical interval count exceeds addressable source data")?;
        Ok(Self { width, source_intervals, output_samples })
    }

    pub(crate) const fn source_intervals(self) -> usize { self.source_intervals }

    pub(crate) fn output_bounds(self, sample: u64) -> Result<(f64, f64), &'static str> {
        if sample >= self.output_samples { return Err("resolved observer sample window is exhausted"); }
        Ok((sample as f64, (sample + 1) as f64))
    }

    pub(crate) fn source_end(self, accepted: u64) -> f64 { accepted as f64 * self.width }

    /// Fraction of one output interval, plus its new overlap cursor. The
    /// caller advances the actual source owner whenever its last interval ends.
    pub(crate) fn overlap(self, cursor: f64, output_end: f64, accepted: u64)
        -> Result<(f64, f64), &'static str>
    {
        if accepted == 0 || accepted > self.source_intervals as u64 {
            return Err("resolved observer has no covering accepted physical interval");
        }
        let source_start = (accepted - 1) as f64 * self.width;
        let next = output_end.min(self.source_end(accepted));
        let weight = next - cursor;
        if source_start > cursor || !weight.is_finite() || weight <= 0.0 || weight > 1.0 {
            return Err("resolved observer overlap is outside its physical time partition");
        }
        Ok((weight, next))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn project(clock: IntervalClock, values: &[f64]) -> Vec<f64> {
        assert_eq!(values.len(), clock.source_intervals());
        let mut accepted = 1u64;
        (0..clock.output_samples).map(|sample| {
            let (mut cursor, end) = clock.output_bounds(sample).unwrap();
            let mut result = 0.0;
            while cursor < end {
                while clock.source_end(accepted) <= cursor { accepted += 1; }
                let (weight, next) = clock.overlap(cursor, end, accepted).unwrap();
                result += weight * values[accepted as usize - 1];
                cursor = next;
            }
            result
        }).collect()
    }

    #[test]
    fn faster_source_intervals_keep_their_exact_integrated_pressure() {
        let clock = IntervalClock::new(0.75, 1, 3).unwrap();
        assert_eq!(clock.source_intervals(), 4);
        let projected = project(clock, &[1.0, 3.0, 5.0, 7.0]);
        assert_eq!(projected, [1.5, 4.0, 6.5]);
        assert_eq!(projected.iter().sum::<f64>(), 0.75 * (1.0 + 3.0 + 5.0 + 7.0));
    }

    #[test]
    fn slower_source_is_retained_across_observer_intervals() {
        let clock = IntervalClock::new(1.5, 1, 3).unwrap();
        assert_eq!(clock.source_intervals(), 2);
        assert_eq!(project(clock, &[2.0, 8.0]), [2.0, 5.0, 8.0]);
    }

    #[test]
    fn final_partial_interval_and_nonbinary_clocks_preserve_constants() {
        let clock = IntervalClock::new(2.5, 1, 3).unwrap();
        assert_eq!(project(clock, &[2.0, 8.0]), [2.0, 2.0, 5.0]);
        for width in [0.1, 1.0 / 3.0, 0.73, 1.37] {
            let clock = IntervalClock::new(width, 1, 1024).unwrap();
            assert!(clock.source_end(clock.source_intervals as u64) >= 1024.0);
            assert!(clock.source_end(clock.source_intervals as u64 - 1) < 1024.0);
            let out = project(clock, &vec![7.0; clock.source_intervals()]);
            assert!(out.iter().all(|value| (*value - 7.0).abs() <= 8.0 * f64::EPSILON * 7.0));
        }
    }

    #[test]
    fn invalid_clocks_or_uncovered_intervals_are_not_silently_repaired() {
        for dt in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::MIN_POSITIVE] {
            assert!(IntervalClock::new(dt, 1, 1).is_err());
        }
        assert!(IntervalClock::new(1.0, 0, 1).is_err());
        assert!(IntervalClock::new(1.0, 1, 0).is_err());
        assert!(IntervalClock::new(1.0, 1, (1 << 53) + 1).is_err());
        let clock = IntervalClock::new(0.75, 1, 3).unwrap();
        assert!(clock.overlap(0.0, 1.0, 0).is_err());
        assert!(clock.overlap(0.0, 1.0, 2).is_err());
        assert!(clock.overlap(0.75, 1.0, 1).is_err());
        assert!(clock.output_bounds(3).is_err());
    }
}
