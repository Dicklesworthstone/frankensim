//! Measurement discipline: warmup, repetition, order statistics.
//!
//! A benchmark without variance bars is folklore (plan §14.1). Every timed
//! kernel reports median plus relative interquartile dispersion so a
//! reader can tell a measurement from a fluke.

/// Order statistics over one kernel's repetition times (seconds).
#[derive(Debug, Clone)]
pub struct Sample {
    /// Median repetition time.
    pub median: f64,
    /// 25th percentile.
    pub p25: f64,
    /// 75th percentile.
    pub p75: f64,
    /// Fastest repetition.
    pub min: f64,
    /// Slowest repetition.
    pub max: f64,
    /// Relative interquartile dispersion: `(p75 − p25) / median`
    /// (0 when the median is 0).
    pub dispersion: f64,
    /// All repetition times, in measurement order (for trend logs).
    pub times: Vec<f64>,
}

/// Time `f` for `reps` repetitions after `warmup` discarded runs.
pub fn time_reps(f: &mut dyn FnMut(), warmup: usize, reps: usize) -> Sample {
    for _ in 0..warmup {
        f();
    }
    let mut times = Vec::with_capacity(reps.max(1));
    for _ in 0..reps.max(1) {
        let start = std::time::Instant::now();
        f();
        times.push(start.elapsed().as_secs_f64());
    }
    sample_from_times(times)
}

/// Order statistics from raw times (meta-tested against hand calculations).
#[must_use]
pub fn sample_from_times(times: Vec<f64>) -> Sample {
    if times.is_empty() {
        // No measurements: a defined zero sample rather than an index/underflow
        // panic on the empty slice.
        return Sample {
            median: 0.0,
            p25: 0.0,
            p75: 0.0,
            min: 0.0,
            max: 0.0,
            dispersion: 0.0,
            times,
        };
    }
    let mut sorted = times.clone();
    sorted.sort_by(f64::total_cmp);
    let q = |p: f64| -> f64 {
        // Nearest-rank on the sorted sample; deterministic tie-breaking.
        let idx = ((sorted.len() as f64 - 1.0) * p).round() as usize;
        sorted[idx.min(sorted.len() - 1)]
    };
    let median = q(0.5);
    let p25 = q(0.25);
    let p75 = q(0.75);
    let dispersion = if median > 0.0 {
        (p75 - p25) / median
    } else {
        0.0
    };
    Sample {
        median,
        p25,
        p75,
        min: sorted[0],
        max: sorted[sorted.len() - 1],
        dispersion,
        times,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn order_statistics_match_hand_calculation() {
        let s = sample_from_times(vec![5.0, 1.0, 3.0, 2.0, 4.0]);
        assert!((s.median - 3.0).abs() < 1e-15);
        assert!((s.p25 - 2.0).abs() < 1e-15);
        assert!((s.p75 - 4.0).abs() < 1e-15);
        assert!((s.min - 1.0).abs() < 1e-15);
        assert!((s.max - 5.0).abs() < 1e-15);
        assert!((s.dispersion - 2.0 / 3.0).abs() < 1e-12);
    }

    #[test]
    fn empty_times_is_defined_not_a_panic() {
        // regression: sorted.len() - 1 underflowed and indexed an empty slice.
        let s = sample_from_times(vec![]);
        assert!((s.median - 0.0).abs() < 1e-15);
        assert!((s.min - 0.0).abs() < 1e-15 && (s.max - 0.0).abs() < 1e-15);
        assert!(s.times.is_empty());
    }

    #[test]
    fn single_rep_is_degenerate_but_defined() {
        let s = sample_from_times(vec![2.5]);
        assert!((s.median - 2.5).abs() < 1e-15);
        assert!((s.dispersion - 0.0).abs() < 1e-15);
    }

    #[test]
    fn time_reps_runs_warmup_and_reps() {
        let mut count = 0u32;
        let s = time_reps(&mut || count += 1, 3, 5);
        assert_eq!(count, 8, "3 warmup + 5 measured");
        assert_eq!(s.times.len(), 5);
        assert!(s.median >= 0.0);
    }
}
