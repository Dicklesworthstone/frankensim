//! Dimensionless local estimates from changes, never absolute state datums.

/// Both arguments are changes from the SAME accepted physical start. Scaling
/// before adding avoids overflowing an otherwise valid absolute-plus-relative
/// denominator. An infinite ratio rejects; it must never become a zero error.
pub(crate) fn scaled_error(coarse: f64, fine: f64, absolute: f64, relative: f64) -> Option<f64> {
    if ![coarse, fine, absolute, relative].into_iter().all(f64::is_finite)
        || absolute <= 0.0 || !(0.0..1.0).contains(&relative)
    {
        return None;
    }
    let change = coarse.abs().max(fine.abs());
    let scale = absolute.max(change);
    let difference = coarse - fine;
    let numerator = if difference.is_finite() {
        difference.abs() / scale
    } else {
        (coarse / scale - fine / scale).abs()
    };
    let denominator = absolute / scale + relative * (change / scale);
    if numerator == 0.0 { Some(0.0) } else { Some(numerator / denominator) }
}

/// The existing split is first order: its coarse/fine local difference is
/// O(dt^2). This only chooses a next interval; no endpoint is extrapolated.
pub(crate) fn step_factor(error: f64, accepted: bool) -> f64 {
    let raw = if error == 0.0 { 2.0 } else { 0.9 / error.sqrt() };
    raw.clamp(0.2, if accepted { 2.0 } else { 0.5 })
}

/// Clamp a suggested comparison duration BEFORE any final-target clipping.
/// Callers admit positive ordered bounds and a positive bounded factor. The
/// division before multiplication also prevents overflow on an allowed growth.
pub(crate) fn bounded_step_s(step_s: f64, factor: f64, minimum_s: f64, maximum_s: f64) -> f64 {
    (step_s.min(maximum_s / factor) * factor).max(minimum_s)
}

/// At a large absolute clock, the two representable half durations can differ.
/// For a first-order method the fine LTE / coarse-fine difference is
/// (a^2+b^2)/(2ab), not invariably one. This form avoids squaring durations.
pub(crate) fn fine_error_factor(first_s: f64, second_s: f64) -> Option<f64> {
    if ![first_s, second_s].into_iter().all(|x| x.is_finite() && x > 0.0) {
        return None;
    }
    let ratio = first_s.max(second_s) / first_s.min(second_s);
    let factor = 0.5 * ratio + 0.5 / ratio;
    factor.is_finite().then_some(factor.max(1.0))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn relative_scale_tracks_changes_and_not_a_large_storage_datum() {
        let datum = 1_000_000_000.0;
        let ratio = scaled_error((datum + 2.0) - datum, (datum + 1.0) - datum, 0.25, 0.25).unwrap();
        assert_eq!(ratio, 1.0 / 0.75);
        assert!(ratio > 1.0);
        assert_eq!(scaled_error(2.0, 1.0, 0.25, 0.25), Some(ratio));
    }

    #[test]
    fn unit_conversion_and_extreme_scales_keep_the_same_admission() {
        assert_eq!(scaled_error(2.0, 1.0, 0.25, 0.25),
            scaled_error(2000.0, 1000.0, 250.0, 0.25));
        let large = scaled_error(f64::MAX, -f64::MAX, f64::MAX, 0.5).unwrap();
        assert!((large - 4.0 / 3.0).abs() < 1e-15);
        assert!(scaled_error(1.0, 0.0, f64::MIN_POSITIVE, 0.0).unwrap() > 1.0);
        assert_eq!(scaled_error(0.0, 0.0, f64::MIN_POSITIVE, 0.0), Some(0.0));
        assert!(scaled_error(f64::NAN, 0.0, 1.0, 0.0).is_none());
        assert!(scaled_error(0.0, 0.0, 0.0, 0.0).is_none());
    }

    #[test]
    fn bounded_first_order_controller_shrinks_every_rejection() {
        assert_eq!(step_factor(0.0, true), 2.0);
        assert_eq!(step_factor(f64::INFINITY, false), 0.2);
        assert!(step_factor(1.01, false) <= 0.5);
        assert_eq!(step_factor(1.0, true), 0.9);
    }

    #[test]
    fn accepted_and_rejected_suggestions_try_the_minimum_without_crossing_it() {
        assert_eq!(bounded_step_s(1.0, step_factor(1.0, true), 1.0, 4.0), 1.0);
        assert_eq!(bounded_step_s(1.5, step_factor(4.0, false), 1.0, 1.5), 1.0);
        assert_eq!(bounded_step_s(1.0, step_factor(4.0, false), 1.0, 1.0), 1.0);
        assert_eq!(bounded_step_s(f64::MAX, 2.0, 1.0, f64::MAX), f64::MAX);
    }

    #[test]
    fn rounded_half_clocks_use_the_actual_fine_truncation_factor() {
        assert_eq!(fine_error_factor(1.0, 1.0), Some(1.0));
        assert_eq!(fine_error_factor(1.0, 2.0), Some(1.25));
        assert_eq!(fine_error_factor(2.0, 1.0), Some(1.25));
        let start = 2.0_f64.powi(52);
        let end = start + 3.0;
        let middle = start + 0.5 * (end - start);
        assert_eq!(fine_error_factor(middle - start, end - middle), Some(1.25));
        assert_eq!(fine_error_factor(f64::MIN_POSITIVE, f64::MIN_POSITIVE), Some(1.0));
        assert_eq!(fine_error_factor(1.0, 0.0), None);
    }
}
