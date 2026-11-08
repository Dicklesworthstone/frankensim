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

/// Successive source-emission roots on one admitted receiver history.
///
/// The original finite source interval remains authoritative. The supplied
/// upper bound on d(tau)/dt only proposes a smaller bracket; both endpoint
/// residuals must confirm it before the existing bisection uses it. A missed
/// prediction, repeated query, or backward query falls back to the full interval.
/// Every solve uses at most `maximum_iterations + 4` residual queries, including
/// the two prediction probes before a fallback. Failed queries retain state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct RetardedEmissionCursor {
    emission_interval_s: [f64; 2],
    maximum_source_time_rate: f64,
    previous: Option<(f64, [f64; 2])>,
}

impl RetardedEmissionCursor {
    /// Parameters are checked by `solve` before querying the original history.
    pub(super) const fn new(emission_interval_s: [f64; 2], maximum_source_time_rate: f64) -> Self {
        Self {
            emission_interval_s,
            maximum_source_time_rate,
            previous: None,
        }
    }

    pub(super) fn solve<E>(
        &mut self,
        reception_s: f64,
        tolerance_s: f64,
        maximum_iterations: u32,
        mut residual_at: impl FnMut(f64) -> Result<f64, E>,
    ) -> Result<[f64; 2], RetardedClockError<E>> {
        let [first, last] = self.emission_interval_s;
        if !(first.is_finite() && last.is_finite() && first <= last)
            || !reception_s.is_finite()
            || !(self.maximum_source_time_rate > 0.0 && self.maximum_source_time_rate.is_finite())
            || !(tolerance_s > 0.0 && tolerance_s.is_finite())
            || !(1..=256).contains(&maximum_iterations)
        {
            return Err(RetardedClockError::InvalidBracket);
        }
        let mut accepted = None;
        if let Some((previous_time, [low, high])) = self.previous
            && reception_s > previous_time
        {
            let elapsed = (reception_s - previous_time).next_up();
            let advance = (elapsed * self.maximum_source_time_rate).next_up();
            let high = (high + advance).next_up().min(last);
            if low <= high && (low > first || high < last) {
                let low_residual = residual_at(low).map_err(RetardedClockError::Query)?;
                let high_residual = residual_at(high).map_err(RetardedClockError::Query)?;
                if !(low_residual.is_finite() && high_residual.is_finite()) {
                    return Err(RetardedClockError::NonFinite);
                }
                if low_residual <= 0.0 && high_residual >= 0.0 {
                    accepted = Some(bisect_root_bracket(
                        low,
                        high,
                        tolerance_s,
                        maximum_iterations,
                        |time| {
                            // Reuse the endpoint probes; neither is an extra
                            // trajectory query or an uncharged retry.
                            if time == low {
                                Ok(low_residual)
                            } else if time == high {
                                Ok(high_residual)
                            } else {
                                residual_at(time)
                            }
                        },
                    )?);
                }
            }
        }
        let bracket = match accepted {
            Some(bracket) => bracket,
            None => bisect_root_bracket(first, last, tolerance_s, maximum_iterations, residual_at)?,
        };
        self.previous = Some((reception_s, bracket));
        Ok(bracket)
    }
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
    if low_residual == 0.0 {
        return Ok([low, low]);
    }
    if high_residual == 0.0 {
        return Ok([high, high]);
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
    fn cursor_preserves_accelerating_retarded_pressure_with_fewer_queries() {
        let mut cursor = RetardedEmissionCursor::new([0.0, 10.0], 1.3);
        let mut cursor_queries = 0_u64;
        let mut full_queries = 0_u64;
        let source_position = |tau: f64| 17.0 * tau + 0.25 * tau * tau;
        for frame in 0..4_096 {
            let reception = 5.0 + f64::from(frame) / 96_000.0;
            let receiver_position = 1_000.0 - 34.0 * reception + reception * reception;
            let residual =
                |tau: f64| tau + (receiver_position - source_position(tau)) / 340.0 - reception;
            let full = bisect_root_bracket(0.0, 10.0, 1.0e-12, 64, |tau| {
                full_queries += 1;
                Ok::<_, ()>(residual(tau))
            })
            .unwrap();
            let before = cursor_queries;
            let narrowed = cursor
                .solve(reception, 1.0e-12, 64, |tau| {
                    cursor_queries += 1;
                    Ok::<_, ()>(residual(tau))
                })
                .unwrap();
            assert!(cursor_queries - before <= 68);
            assert!(narrowed[0] >= 0.0 && narrowed[1] <= 10.0);
            assert!(residual(narrowed[0]) <= 0.0);
            assert!(residual(narrowed[1]) >= 0.0);
            let emission = 0.5 * (narrowed[0] + narrowed[1]);
            let reference = 0.5 * (full[0] + full[1]);
            // Independent smaller quadratic root for these accelerating paths.
            let q = 374.0 * reception - 1_000.0 - reception * reception;
            let exact = 2.0 * q / (323.0 + (323.0 * 323.0 - q).sqrt());
            assert!((emission - exact).abs() <= 1.0e-12);
            assert!((emission - reference).abs() <= 1.0e-12);
            let pressure = |tau: f64| {
                (core::f64::consts::TAU * 1_200.0 * tau).sin()
                    / (receiver_position - source_position(tau))
            };
            assert!((pressure(emission) - pressure(reference)).abs() < 1.0e-10);
        }
        assert!(cursor_queries * 4 < full_queries * 3);
        eprintln!("4096 physical clock queries: cursor={cursor_queries}, full={full_queries}");
    }

    #[test]
    fn cursor_checks_predictions_and_preserves_state_on_failed_queries() {
        // Deliberately understate the clock-rate hint. It must cause a verified
        // full-domain fallback, never clipping or a false successful root.
        let mut cursor = RetardedEmissionCursor::new([0.0, 10.0], 0.25);
        for reception in [1.0, 2.0, 2.0, 1.5, 1.500_1] {
            let expected = 0.75 + 2.0 * reception;
            let mut queries = 0;
            let bracket = cursor
                .solve(reception, 1.0e-12, 64, |tau| {
                    queries += 1;
                    Ok::<_, &str>(tau - expected)
                })
                .unwrap();
            assert!(bracket[0] <= expected && expected <= bracket[1]);
            assert!(bracket[1] - bracket[0] <= 1.0e-12);
            assert!(queries <= 68);
        }
        let saved = cursor;
        let mut queries = 0;
        assert_eq!(
            cursor.solve(1.51, 1.0e-12, 64, |tau| {
                queries += 1;
                if queries == 3 {
                    Err("cancelled")
                } else {
                    Ok(tau - 3.77)
                }
            }),
            Err(RetardedClockError::Query("cancelled"))
        );
        assert_eq!(queries, 3);
        assert_eq!(cursor, saved);
        assert_eq!(
            cursor.solve(8.0, 1.0e-12, 64, |tau| Ok::<_, ()>(tau - 16.75)),
            Err(RetardedClockError::InvalidBracket)
        );
        assert_eq!(cursor, saved);
        assert_eq!(
            cursor.solve(1.51, 1.0e-12, 64, |_| Ok::<_, ()>(f64::NAN)),
            Err(RetardedClockError::NonFinite)
        );
        assert_eq!(cursor, saved);
        assert_eq!(
            cursor.solve(1.51, 1.0e-12, 1, |tau| Ok::<_, ()>(tau - 3.77)),
            Err(RetardedClockError::IterationLimit)
        );
        assert_eq!(cursor, saved);
    }

    #[test]
    fn exact_endpoint_roots_do_not_consume_the_iteration_budget() {
        for root in [0.0, 1.0] {
            let mut queries = 0;
            let bracket = bisect_root_bracket(0.0, 1.0, 1.0e-12, 1, |time| {
                queries += 1;
                Ok::<_, ()>(time - root)
            })
            .unwrap();
            assert_eq!(bracket, [root, root]);
            assert_eq!(queries, 2);
        }
        for rate in [0.0, f64::NAN, f64::INFINITY] {
            let mut cursor = RetardedEmissionCursor::new([0.0, 1.0], rate);
            let mut queries = 0;
            assert_eq!(
                cursor.solve(0.5, 1.0e-12, 64, |time| {
                    queries += 1;
                    Ok::<_, ()>(time - 0.5)
                }),
                Err(RetardedClockError::InvalidBracket)
            );
            assert_eq!(queries, 0);
        }
    }

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
                })
                .expect("subsonic emission");
                let emission = 0.5 * (low + high);
                let expected =
                    (arrival * (1.0 - receiver_v / c) - separation / c) / (1.0 - source_v / c);
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
        let residual =
            |endpoint: usize, time: f64| time - (emission[endpoint] + position(time) / c);
        let window =
            clip_reception_window([0.0, 2.0], 1.0e-12, 64, |e, t| Ok::<_, ()>(residual(e, t)))
                .expect("finite waves");
        for (endpoint, boundary) in window.into_iter().enumerate() {
            let expected = (330.0_f64
                - (330.0_f64 * 330.0 - 8.0 * (100.0 + c * emission[endpoint])).sqrt())
                / 4.0;
            assert!((boundary - expected).abs() < 2.0e-12);
        }
        assert!(residual(0, window[0]) >= 0.0);
        assert!(residual(1, window[1]) <= 0.0);
        assert!((window[0] - position(0.0) / c).abs() > 0.005);
        let clipped =
            clip_reception_window([0.5, 1.0], 1.0e-12, 64, |e, t| Ok::<_, ()>(residual(e, t)))
                .expect("receiver subset is retained exactly");
        assert_eq!(clipped, [0.5, 1.0]);
    }

    #[test]
    fn finite_histories_refuse_disjoint_windows_and_queries() {
        for interval in [[0.0, 0.5], [3.0, 4.0]] {
            assert!(
                clip_reception_window(interval, 1.0e-12, 64, |endpoint, time| {
                    Ok::<_, ()>(time - (1.0 + endpoint as f64))
                })
                .is_err()
            );
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
