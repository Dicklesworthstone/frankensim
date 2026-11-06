//! Bounded scalar clocks shared by source emission and receiver-history cuts.
//!
//! Motion admission belongs to the original trajectory owners. Their strictly
//! subsonic bounds make these residuals monotone; this module only inverts them.

#[derive(Debug, PartialEq)]
pub(super) enum RetardedClockError<E> {
    Query(E),
    InvalidBracket,
    NonFinite,
    IterationLimit,
}

/// Return an enclosing bracket; an exactly represented root collapses it.
/// Callers choosing a finite-history cut select the appropriate INNER endpoint.
pub(super) fn bisect_root_bracket<E>(
    mut low: f64,
    mut high: f64,
    tolerance_s: f64,
    maximum_iterations: u32,
    mut residual_at: impl FnMut(f64) -> Result<f64, E>,
) -> Result<[f64; 2], RetardedClockError<E>> {
    if !(low.is_finite() && high.is_finite() && low <= high)
        || !(tolerance_s > 0.0 && tolerance_s.is_finite())
        || !(1..=256).contains(&maximum_iterations)
    {
        return Err(RetardedClockError::InvalidBracket);
    }
    let mut low_residual = residual_at(low).map_err(RetardedClockError::Query)?;
    let high_residual = residual_at(high).map_err(RetardedClockError::Query)?;
    if !(low_residual.is_finite() && high_residual.is_finite()) {
        return Err(RetardedClockError::NonFinite);
    }
    if low_residual > 0.0 || high_residual < 0.0 {
        return Err(RetardedClockError::InvalidBracket);
    }
    for _ in 0..maximum_iterations {
        let mid = 0.5 * (low + high);
        let residual = residual_at(mid).map_err(RetardedClockError::Query)?;
        if !residual.is_finite() {
            return Err(RetardedClockError::NonFinite);
        }
        if residual <= 0.0 {
            low = mid;
            low_residual = residual;
        } else {
            high = mid;
        }
        if high - low <= tolerance_s || residual == 0.0 {
            break;
        }
    }
    if high - low > tolerance_s && low_residual != 0.0 {
        return Err(RetardedClockError::IterationLimit);
    }
    Ok(if low_residual == 0.0 {
        [low, low]
    } else {
        [low, high]
    })
}

/// Intersect finite receiver history with the waves from the first and last
/// complete source stencils. `residual_at(endpoint, t)` is
/// `t - (tau_endpoint + distance(receiver(t), source(tau_endpoint))/c)`.
/// It increases strictly for a subsonic receiver. The first wave must already
/// have arrived; the last must not yet have arrived. No history is extrapolated.
pub(super) fn clip_reception_window<E>(
    receiver_interval_s: [f64; 2],
    tolerance_s: f64,
    maximum_iterations: u32,
    mut residual_at: impl FnMut(usize, f64) -> Result<f64, E>,
) -> Result<[f64; 2], RetardedClockError<E>> {
    let [first, last] = receiver_interval_s;
    if !(first.is_finite() && last.is_finite() && first < last) {
        return Err(RetardedClockError::InvalidBracket);
    }
    let first_wave_at_start = residual_at(0, first).map_err(RetardedClockError::Query)?;
    let last_wave_at_end = residual_at(1, last).map_err(RetardedClockError::Query)?;
    if !(first_wave_at_start.is_finite() && last_wave_at_end.is_finite()) {
        return Err(RetardedClockError::NonFinite);
    }
    let lower = if first_wave_at_start >= 0.0 {
        first
    } else {
        bisect_root_bracket(first, last, tolerance_s, maximum_iterations, |time| {
            residual_at(0, time)
        })?[1]
    };
    let upper = if last_wave_at_end <= 0.0 {
        last
    } else {
        bisect_root_bracket(first, last, tolerance_s, maximum_iterations, |time| {
            residual_at(1, time)
        })?[0]
    };
    if lower >= upper {
        return Err(RetardedClockError::InvalidBracket);
    }
    Ok([lower, upper])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn both_endpoint_velocities_set_the_retarded_clock_rate() {
        let c = 340.0;
        let separation = 100.0;
        for (source_v, receiver_v) in [(17.0, -34.0), (-17.0, 34.0), (17.0, 17.0)] {
            let residual = |tau: f64, t: f64| {
                let radius = separation + receiver_v * t - source_v * tau;
                assert!(radius > 0.0);
                tau + radius / c - t
            };
            let expected_rate = (1.0 - receiver_v / c) / (1.0 - source_v / c);
            let mut times = Vec::new();
            for arrival in [1.0, 1.01] {
                let [low, high] = bisect_root_bracket(0.0, 2.0, 1.0e-13, 64, |tau| {
                    Ok::<_, ()>(residual(tau, arrival))
                }).expect("subsonic emission");
                let emission = 0.5 * (low + high);
                let expected = (arrival * (1.0 - receiver_v / c) - separation / c)
                    / (1.0 - source_v / c);
                assert!((emission - expected).abs() < 1.0e-12);
                times.push(emission);
            }
            assert!(((times[1] - times[0]) / 0.01 - expected_rate).abs() < 2.0e-11);
        }
    }

    #[test]
    fn accelerating_receiver_cuts_use_wavefront_roots_and_inner_brackets() {
        let c = 340.0;
        let emission = [0.0, 1.0];
        let position = |time: f64| 100.0 + 10.0 * time + 2.0 * time * time;
        let residual = |endpoint: usize, time: f64| {
            time - (emission[endpoint] + position(time) / c)
        };
        let window = clip_reception_window([0.0, 2.0], 1.0e-12, 64, |e, t| {
            Ok::<_, ()>(residual(e, t))
        }).expect("finite waves");
        for (endpoint, boundary) in window.into_iter().enumerate() {
            let expected = (330.0_f64
                - (330.0_f64 * 330.0 - 8.0 * (100.0 + c * emission[endpoint])).sqrt())
                / 4.0;
            assert!((boundary - expected).abs() < 2.0e-12);
        }
        assert!(residual(0, window[0]) >= 0.0);
        assert!(residual(1, window[1]) <= 0.0);
        assert!((window[0] - position(0.0) / c).abs() > 0.005);
        let clipped = clip_reception_window([0.5, 1.0], 1.0e-12, 64, |e, t| {
            Ok::<_, ()>(residual(e, t))
        }).expect("receiver subset is retained exactly");
        assert_eq!(clipped, [0.5, 1.0]);
    }

    #[test]
    fn finite_histories_refuse_disjoint_windows_and_queries() {
        for interval in [[0.0, 0.5], [3.0, 4.0]] {
            assert!(clip_reception_window(interval, 1.0e-12, 64, |endpoint, time| {
                Ok::<_, ()>(time - (1.0 + endpoint as f64))
            }).is_err());
        }
        assert_eq!(
            clip_reception_window([0.0, 4.0], 1.0e-12, 64, |_, _| Err("cancelled")),
            Err(RetardedClockError::Query("cancelled"))
        );
        assert_eq!(
            bisect_root_bracket(0.0, 1.0, 1.0e-15, 1, |x| Ok::<_, ()>(x - 0.3)),
            Err(RetardedClockError::IterationLimit)
        );
        assert_eq!(
            bisect_root_bracket(0.0, 1.0, 1.0e-12, 64, |_| Ok::<_, ()>(f64::NAN)),
            Err(RetardedClockError::NonFinite)
        );
    }
}
