//! Positive-feedback bound for one original affine electrical-source interval.
//! Native conductor powers are supplied by the graph owner, never re-solved here.

/// For q(T)>0, P=C*q has constant derivative; P=C/q has its largest
/// positive derivative at the smallest q. Only the declared source is affine,
/// not its reciprocal. Inputs are source endpoints and their native powers.
pub(super) fn positive_segment_slope_w_k(
    width_k: f64,
    source: [f64; 2],
    power_w: [f64; 2],
    inverse: bool,
) -> Option<f64> {
    if !width_k.is_finite()
        || width_k <= 0.0
        || source.iter().any(|q| !q.is_finite() || *q <= 0.0)
        || power_w.iter().any(|p| !p.is_finite() || *p < 0.0)
    {
        return None;
    }
    let change = source[1] - source[0];
    if change == 0.0 || (inverse && change > 0.0) || (!inverse && change < 0.0) {
        return Some(0.0);
    }
    // Positive feedback always attains its maximum at the right endpoint:
    // q increases for the direct law and decreases for the reciprocal law.
    let power = power_w[1];
    if power == 0.0 {
        return Some(0.0);
    }
    let delta = change.abs().next_up();
    let width = width_k.next_down();
    if width <= 0.0 || !delta.is_finite() {
        return None;
    }
    let q = source[1];
    let multiply = |a: f64, b: f64| (a * b).next_up();
    let divide = |a: f64, b: f64| (a / b).next_up();
    // Equivalent positive factorizations avoid rejecting a finite bound just
    // because P/q or a numerator product exceeds the intermediate range.
    let bound = [
        divide(multiply(divide(power, q), delta), width),
        divide(multiply(power, divide(delta, q)), width),
        multiply(divide(power, width), divide(delta, q)),
        divide(multiply(power, divide(delta, width)), q),
    ]
    .into_iter()
    .filter(|value| value.is_finite())
    .min_by(f64::total_cmp)?;
    // Retain a roundoff allowance for the canonical source, resistance,
    // current and power evaluations that precede this interval algebra.
    let bound = (bound * (1.0 + 64.0 * f64::EPSILON)).next_up();
    bound.is_finite().then_some(bound)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bounded(actual: f64, exact: f64) {
        assert!(actual >= exact, "{actual} < {exact}");
        assert!(actual <= exact * (1.0 + 1.0e-12), "{actual} > {exact}");
    }

    #[test]
    fn g1_direct_and_reciprocal_bounds_follow_the_original_source_law() {
        bounded(
            positive_segment_slope_w_k(100.0, [2.0, 3.0], [72.0, 108.0], false).unwrap(),
            0.36,
        );
        let reciprocal = positive_segment_slope_w_k(100.0, [3.0, 2.0], [48.0, 72.0], true).unwrap();
        bounded(reciprocal, 0.36);
        // Independent analytic derivative at interior temperatures. Using
        // a secant of reciprocal endpoint powers would underestimate this.
        for temperature in [300.0, 325.0, 350.0, 375.0, 400.0] {
            let q = 3.0 - 0.01 * (temperature - 300.0);
            assert!(144.0 * 0.01 / (q * q) <= reciprocal);
        }
    }

    #[test]
    fn g0_nonpositive_feedback_and_zero_drive_have_zero_bound() {
        for (source, power, inverse) in [
            ([2.0, 3.0], [72.0, 48.0], true),
            ([3.0, 2.0], [108.0, 72.0], false),
            ([2.0, 2.0], [72.0, 72.0], false),
            ([2.0, 3.0], [0.0, 0.0], false),
        ] {
            assert_eq!(
                positive_segment_slope_w_k(100.0, source, power, inverse),
                Some(0.0)
            );
        }
    }

    #[test]
    fn g3_source_scaling_preserves_the_same_power_response() {
        for scale in [1.0e-300, 1.0, 1.0e300] {
            bounded(
                positive_segment_slope_w_k(100.0, [3.0 * scale, 2.0 * scale], [48.0, 72.0], true)
                    .unwrap(),
                0.36,
            );
        }
    }

    #[test]
    fn g0_tiny_positive_feedback_never_rounds_to_zero() {
        let bound = positive_segment_slope_w_k(
            100.0,
            [1.0, 1.0 + f64::EPSILON],
            [1.0e-300, 1.0e-300 * (1.0 + f64::EPSILON)],
            false,
        )
        .unwrap();
        assert!(bound > 0.0 && bound.is_finite());
    }

    #[test]
    fn g0_invalid_intervals_or_power_do_not_become_zero_feedback() {
        for (width, source, power) in [
            (0.0, [2.0, 3.0], [1.0, 1.5]),
            (100.0, [0.0, 3.0], [1.0, 1.5]),
            (100.0, [2.0, 2.0], [1.0, f64::INFINITY]),
            (100.0, [2.0, 3.0], [1.0, f64::NAN]),
        ] {
            assert_eq!(
                positive_segment_slope_w_k(width, source, power, false),
                None
            );
        }
    }
}
