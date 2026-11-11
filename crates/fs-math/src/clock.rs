//! Rounded simulation clocks that must land exactly on their end time.
//!
//! A scheduler that advances `elapsed = elapsed + step` toward `target` rounds
//! once per addition. Steps whose exact sum is `target - start` can then stop a
//! few float spacings short of `target` (ten 0.05 s steps end at
//! 0.49999999999999994, not 0.5). The interval left behind is too short to
//! carry any resolvable state change, and the next step on it refuses.
//! [`clock_endpoint`] runs the step that would leave such an interval on to
//! `target` instead. The allowance is the clock's own rounding bound, not a
//! tolerance on the physics.

use crate::canonical_nan;

/// Unit roundoff `u = 2^-53` of binary64 round-to-nearest.
const UNIT_ROUNDOFF: f64 = f64::EPSILON / 2.0;

/// Endpoint of one clock step from `elapsed` toward `target`.
///
/// `start` is the clock's value at the beginning of the caller's interval, and
/// `additions` is the number of rounded additions that produced `elapsed` from
/// it (0 for the first step). The result is `target` when the step reaches or
/// passes it, or when the interval the step would leave,
/// `target - (elapsed + step)`, is at most `(additions + 2) u M` with
/// `M = max(|start|, |target|)`. Otherwise it is the rounded `elapsed + step`.
///
/// The bound. Each round-to-nearest addition gives `fl(a + b) = (a + b) /
/// (1 + d)` with `|d| <= u` (Higham, Accuracy and Stability of Numerical
/// Algorithms, 2nd ed., §2.2), so it errs by at most `u |fl(a + b)|`. Summed
/// recursively, the computed clock differs from `start` plus the exact sum of
/// its steps by at most `u` times the sum of its partial-sum magnitudes (§4.2).
/// Every partial sum lies between `start` and `target`, so after this step's
/// endpoint the clock is within `(additions + 1) u M` of where exact steps would
/// put it. A leftover of one float spacing below `target` cannot be divided by
/// any later step either; it is at most `ulp(target) <= 2 u |target| <= 2 u M`.
/// `(additions + 2) u M` covers both cases. With `additions = 0` the allowance
/// is `2 u M`: the one-spacing rule alone.
///
/// The final step then exceeds `step` by at most `(additions + 2) u M`, a
/// clock-rounding amount that callers state beside their step bounds. A leftover
/// above the allowance is resolvable and is left for the caller's next step.
///
/// Requires finite `start <= elapsed < target` and `step > 0`; otherwise the
/// canonical NaN is returned. The count is converted to f64; above 2^53
/// additions it rounds by a relative 2^-53, far below the allowance's meaning.
#[must_use]
pub fn clock_endpoint(start: f64, elapsed: f64, step: f64, target: f64, additions: usize) -> f64 {
    let ordered = start.is_finite()
        && elapsed.is_finite()
        && target.is_finite()
        && start <= elapsed
        && elapsed < target
        && step > 0.0;
    if !ordered {
        return canonical_nan();
    }
    if step >= target - elapsed {
        return target;
    }
    let endpoint = elapsed + step;
    #[allow(clippy::cast_precision_loss)] // see the conversion note above
    let count = additions as f64;
    let allowance = (count + 2.0) * UNIT_ROUNDOFF * start.abs().max(target.abs());
    if endpoint >= target || target - endpoint <= allowance {
        target
    } else {
        endpoint
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const U: f64 = UNIT_ROUNDOFF;

    /// The clock after `count` additions of `step` from `start`.
    fn summed(start: f64, step: f64, count: usize) -> f64 {
        let mut clock = start;
        for _ in 0..count {
            clock += step;
        }
        clock
    }

    #[test]
    fn ten_crystallization_trials_of_a_twentieth_reach_half_a_second() {
        // fs-thermal crystallization step doubling (bead rc-1109 .1.5, B2):
        // the tenth 0.05 s trial ended at 0.49999999999999994, one spacing
        // and 1.0 u t_end short; the eleventh trial's midpoint was its end.
        let clock = summed(0.0, 0.05, 9);
        assert_eq!(clock, 0.449_999_999_999_999_96);
        assert_eq!(clock + 0.05, 0.499_999_999_999_999_94);
        assert_eq!((0.5 - (clock + 0.05)) / (U * 0.5), 1.0);
        assert_eq!(clock_endpoint(0.0, clock, 0.05, 0.5, 9), 0.5);
    }

    #[test]
    fn five_hundredths_after_an_exact_birth_event_reach_the_tenth() {
        // fs-material sublimation timeline (S4, late birth): the clock lands
        // exactly on the 0.05 s birth event, then five 0.01 s steps end at
        // 0.09999999999999999, 1.25 u t_end short of 0.1.
        let clock = summed(0.05, 0.01, 4);
        assert_eq!(clock, 0.09);
        assert_eq!(clock + 0.01, 0.099_999_999_999_999_99);
        assert_eq!((0.1 - (clock + 0.01)) / (U * 0.1), 1.25);
        assert_eq!(clock_endpoint(0.0, clock, 0.01, 0.1, 5), 0.1);
    }

    #[test]
    fn two_thousand_ten_thousandths_reach_a_fifth_of_a_second() {
        // fs-material sublimation timeline (S4, adaptive plateau): 2000 steps
        // of 1e-4 s end at 0.1999999999999943, 257.5 u t_end short. No small
        // fixed multiple of u t_end covers it; the addition count does.
        let clock = summed(0.0, 1e-4, 1999);
        assert_eq!(clock + 1e-4, 0.199_999_999_999_994_3);
        assert_eq!((0.2 - (clock + 1e-4)) / (U * 0.2), 257.5);
        assert!(0.2 - (clock + 1e-4) <= 2001.0 * U * 0.2);
        assert_eq!(clock_endpoint(0.0, clock, 1e-4, 0.2, 1999), 0.2);
    }

    #[test]
    fn a_three_spacing_aeration_leftover_runs_on_only_as_far_as_the_count_allows() {
        // fs-material cyanide aeration (O2): the 16th accepted step of a call
        // from 0.4 s to 1.0 s ended at 0.9999999999999997, three spacings
        // (3.0 u) short. Its midpoint is strictly inside, so a divisibility
        // rule keeps it, but no stage can resolve its Au progress.
        let (elapsed, step) = (0.9, 0.099_999_999_999_999_64);
        assert_eq!(elapsed + step, 0.999_999_999_999_999_7);
        assert_eq!(1.0 - (elapsed + step), 3.0 * U);
        assert_eq!(clock_endpoint(0.4, elapsed, step, 1.0, 15), 1.0);
        // Planted negative: with no prior additions the allowance is 2 u M,
        // the leftover is resolvable by that bound, and the step keeps its own
        // endpoint. One prior addition makes the allowance exactly 3 u M.
        assert_eq!(
            clock_endpoint(0.4, elapsed, step, 1.0, 0),
            0.999_999_999_999_999_7
        );
        assert_eq!(clock_endpoint(0.4, elapsed, step, 1.0, 1), 1.0);
    }

    #[test]
    fn no_prior_additions_run_on_only_one_spacing() {
        // k = 0 is the one-spacing rule: a step ending one spacing below the
        // target runs on, and one ending four spacings below keeps its end.
        let one = 0.5_f64.next_down();
        assert_eq!(clock_endpoint(0.0, 0.0, one, 0.5, 0), 0.5);
        let four = 1.0 - 4.0 * U;
        assert_eq!(clock_endpoint(0.0, 0.0, four, 1.0, 0), four);
        assert_eq!(clock_endpoint(0.0, 0.0, 0.05, 0.5, 0), 0.05);
    }

    #[test]
    fn a_clock_from_below_zero_is_bounded_by_its_start() {
        // 23 equal steps from -2.0 to 0.1 end 52.5 u |target| short of 0.1,
        // above (22 + 2) u |target| but within (22 + 2) u max(|start|, |target|).
        let step = (0.1 - -2.0) / 23.0;
        let clock = summed(-2.0, step, 22);
        let left = 0.1 - (clock + step);
        assert_eq!(left / (U * 0.1), 52.5);
        assert!(left > 24.0 * U * 0.1 && left <= 24.0 * U * 2.0);
        assert_eq!(clock_endpoint(-2.0, clock, step, 0.1, 22), 0.1);
    }

    #[test]
    fn a_step_that_reaches_the_target_ends_on_it_and_bad_input_is_nan() {
        assert_eq!(clock_endpoint(0.0, 0.25, 0.25, 0.5, 3), 0.5);
        assert_eq!(clock_endpoint(0.0, 0.25, 1.0, 0.5, 3), 0.5);
        assert_eq!(clock_endpoint(0.0, 0.25, f64::INFINITY, 0.5, 3), 0.5);
        for (start, elapsed, step, target) in [
            (0.0, 0.5, 0.1, 0.5),
            (0.0, 0.6, 0.1, 0.5),
            (0.3, 0.2, 0.1, 0.5),
            (0.0, 0.1, 0.0, 0.5),
            (0.0, 0.1, -0.1, 0.5),
            (0.0, 0.1, f64::NAN, 0.5),
            (f64::NEG_INFINITY, 0.1, 0.1, 0.5),
            (0.0, 0.1, 0.1, f64::INFINITY),
        ] {
            assert!(clock_endpoint(start, elapsed, step, target, 0).is_nan());
        }
    }
}
