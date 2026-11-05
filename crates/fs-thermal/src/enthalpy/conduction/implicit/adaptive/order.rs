//! Method-order arithmetic used by the native step-doubling controller.
//!
//! This module has no solver or constitutive dependencies. Both supported
//! methods accept the conservative fine trajectory, never an extrapolated
//! state. The correction estimates its local error from the coarse/fine
//! difference using the actual two representable fine spans.

/// Temporal order declared by a native endpoint adapter.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ConductionTimeOrder {
    /// Backward Euler and existing adapters retain their original controller.
    #[default]
    First,
    /// Native SDIRK2 on a smooth admitted solution; corners may reduce order.
    Second,
}

impl ConductionTimeOrder {
    /// Fine local error divided by the coarse/fine discrepancy. For spans a,b,
    /// this is (a²+b²)/(2ab) at order one and
    /// (a³+b³)/(3ab(a+b)) at order two. Equal halves give 1 and 1/3.
    ///
    /// The second-order expression uses the equivalent dimensionless ratio
    /// (1-r+r²)/(3r), r=min(a,b)/max(a,b), without cubing physical seconds.
    /// Invalid spans or an unrepresentable positive correction return None.
    pub fn fine_error_correction(self, a: f64, b: f64) -> Option<f64> {
        if !a.is_finite() || !b.is_finite() || a <= 0.0 || b <= 0.0 {
            return None;
        }
        let correction = match self {
            // Preserve the original first-order arithmetic exactly.
            Self::First => 0.5 * (a / b + b / a),
            Self::Second => {
                let ratio = a.min(b) / a.max(b);
                (1.0 - ratio + ratio * ratio) / (3.0 * ratio)
            }
        };
        (correction.is_finite() && correction > 0.0).then_some(correction)
    }

    /// Suggested multiplier after accepting a normalized local error.
    /// The original safety factor 0.9 and growth bounds [0.5,2] are retained.
    /// Positive infinity requests the lower bound; NaN/negative errors refuse.
    pub fn accepted_growth_factor(self, normalized_error: f64) -> Option<f64> {
        self.controller_factor(normalized_error, 0.5, 2.0)
    }

    /// Suggested multiplier after an accuracy rejection, bounded to [0.1,0.5].
    /// Positive infinity requests the strongest allowed reduction.
    pub fn rejected_reduction_factor(self, normalized_error: f64) -> Option<f64> {
        self.controller_factor(normalized_error, 0.1, 0.5)
    }

    fn controller_factor(self, error: f64, lower: f64, upper: f64) -> Option<f64> {
        if error.is_nan() || error < 0.0 {
            return None;
        }
        Some(match self {
            Self::First => {
                // Keep both the zero-error case and IEEE sqrt path unchanged.
                let factor = if error == 0.0 {
                    2.0
                } else {
                    0.9 / error.sqrt()
                };
                factor.clamp(lower, upper)
            }
            Self::Second => cubic_factor(error, lower, upper),
        })
    }
}

// Solve error*factor³ = 0.9³ inside the already bounded controller interval.
// Comparing error with (0.9/factor)³ keeps every intermediate finite even
// for an infinite normalized error. Fixed IEEE operations and bounded bisection
// avoid a platform cbrt/powf in deterministic timestep selection.
fn cubic_factor(error: f64, mut lower: f64, mut upper: f64) -> f64 {
    let threshold = |factor: f64| {
        let ratio = 0.9 / factor;
        ratio * ratio * ratio
    };
    if error >= threshold(lower) {
        return lower;
    }
    if error <= threshold(upper) {
        return upper;
    }
    for _ in 0..64 {
        let middle = 0.5 * (lower + upper);
        if middle == lower || middle == upper {
            break;
        }
        if error > threshold(middle) {
            upper = middle;
        } else {
            lower = middle;
        }
    }
    lower
}

#[cfg(test)]
mod tests {
    use super::*;

    fn near(actual: f64, expected: f64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance,
            "{actual} != {expected}"
        );
    }

    #[test]
    fn first_order_defaults_and_arithmetic_remain_unchanged() {
        assert_eq!(ConductionTimeOrder::default(), ConductionTimeOrder::First);
        for (a, b) in [(1.0, 1.0), (2.0, 1.0), (0.03125, 0.09375), (1e100, 3e100)] {
            assert_eq!(
                ConductionTimeOrder::First
                    .fine_error_correction(a, b)
                    .unwrap()
                    .to_bits(),
                (0.5 * (a / b + b / a)).to_bits(),
            );
        }
        for error in [0.0, 1e-100, 0.01, 0.3, 1.0, 4.0, 20.0, f64::INFINITY] {
            let raw = if error == 0.0 {
                2.0
            } else {
                0.9 / error.sqrt()
            };
            assert_eq!(
                ConductionTimeOrder::First
                    .accepted_growth_factor(error)
                    .unwrap()
                    .to_bits(),
                raw.clamp(0.5, 2.0).to_bits()
            );
            assert_eq!(
                ConductionTimeOrder::First
                    .rejected_reduction_factor(error)
                    .unwrap()
                    .to_bits(),
                raw.clamp(0.1, 0.5).to_bits()
            );
        }
    }

    #[test]
    fn second_order_uses_actual_unequal_spans_without_dimensional_overflow() {
        let order = ConductionTimeOrder::Second;
        near(
            order.fine_error_correction(1.0, 1.0).unwrap(),
            1.0 / 3.0,
            0.0,
        );
        near(order.fine_error_correction(2.0, 1.0).unwrap(), 0.5, 0.0);
        for (a, b) in [(2.0, 3.0), (1.0, 4.0), (7.0, 9.0)] {
            let expected = (a * a * a + b * b * b) / (3.0 * a * b * (a + b));
            near(order.fine_error_correction(a, b).unwrap(), expected, 4e-16);
            assert_eq!(
                order.fine_error_correction(a, b),
                order.fine_error_correction(b, a)
            );
        }
        for scale in [f64::from_bits(1), 1e-300, 1e300, f64::MAX] {
            near(
                order.fine_error_correction(scale, scale).unwrap(),
                1.0 / 3.0,
                0.0,
            );
        }
        near(order.fine_error_correction(1e300, 5e299).unwrap(), 0.5, 0.0);
        near(
            order.fine_error_correction(1e-300, 5e-301).unwrap(),
            0.5,
            0.0,
        );
    }

    #[test]
    fn cubic_controller_uses_order_three_error_scaling_and_fixed_bounds() {
        let order = ConductionTimeOrder::Second;
        for factor in [0.51, 0.7, 0.9, 1.2, 1.75, 1.99] {
            let ratio = 0.9 / factor;
            let error = ratio * ratio * ratio;
            let actual = order.accepted_growth_factor(error).unwrap();
            near(actual, factor, 8e-16);
            assert_eq!(
                actual.to_bits(),
                order.accepted_growth_factor(error).unwrap().to_bits()
            );
        }
        for factor in [0.11, 0.2, 0.35, 0.49] {
            let ratio = 0.9 / factor;
            near(
                order
                    .rejected_reduction_factor(ratio * ratio * ratio)
                    .unwrap(),
                factor,
                2e-16,
            );
        }
        near(order.rejected_reduction_factor(64.0).unwrap(), 0.225, 2e-16);
        assert_eq!(order.accepted_growth_factor(0.0), Some(2.0));
        assert_eq!(order.accepted_growth_factor(f64::INFINITY), Some(0.5));
        assert_eq!(order.rejected_reduction_factor(f64::INFINITY), Some(0.1));
        assert_eq!(order.rejected_reduction_factor(1.0), Some(0.5));
    }

    #[test]
    fn invalid_spans_and_errors_refuse_without_nan_controller_state() {
        for order in [ConductionTimeOrder::First, ConductionTimeOrder::Second] {
            for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
                assert_eq!(order.fine_error_correction(bad, 1.0), None);
                assert_eq!(order.fine_error_correction(1.0, bad), None);
            }
            for bad in [-1.0, f64::NEG_INFINITY, f64::NAN] {
                assert_eq!(order.accepted_growth_factor(bad), None);
                assert_eq!(order.rejected_reduction_factor(bad), None);
            }
            assert_eq!(
                order.fine_error_correction(f64::MAX, f64::from_bits(1)),
                None
            );
        }
    }
}
