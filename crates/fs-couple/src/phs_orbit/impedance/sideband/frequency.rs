//! One real-source frequency-symmetry seam for cached ports and physical outputs.
use fs_math::c64::C64;

/// `positive` is an exp(+iwt) value or its derivative at |w|. A real source obeys
/// F(-w)=conj(F(w)); differentiating the absolute-frequency argument adds a minus
/// sign. At DC admission already requires a real value and purely imaginary
/// derivative, so both signed limits agree. No transpose or peak scaling occurs.
#[must_use]
pub(crate) fn signed_frequency_entry(positive: C64, frequency: f64, derivative: bool) -> C64 {
    if frequency < 0.0 {
        let reflected = positive.conj();
        if derivative { -reflected } else { reflected }
    } else {
        positive
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn derivative_reflection_changes_the_real_part_and_preserves_the_imaginary_part() {
        let value = C64::new(-0.4, 0.2);
        assert_eq!(
            signed_frequency_entry(value, -2.0, false),
            C64::new(-0.4, -0.2)
        );
        assert_eq!(
            signed_frequency_entry(value, -2.0, true),
            C64::new(0.4, 0.2)
        );
        assert_eq!(signed_frequency_entry(value, 2.0, true), value);
    }

    #[test]
    fn dispersive_original_values_and_slopes_agree_at_both_signed_frequencies() {
        for w in [-11.0_f64, -2.0, -0.3, 0.0, 0.3, 2.0, 11.0] {
            let positive_pole = C64::new(1.0, 0.4 * w.abs());
            let value = C64::new(0.3, 0.2 * w.abs()) + C64::from_re(0.7) / positive_pole;
            let slope = C64::new(0.0, 0.2) + C64::new(0.0, -0.28) / (positive_pole * positive_pole);
            let signed_pole = C64::new(1.0, 0.4 * w);
            let expected = C64::new(0.3, 0.2 * w) + C64::from_re(0.7) / signed_pole;
            let expected_slope =
                C64::new(0.0, 0.2) + C64::new(0.0, -0.28) / (signed_pole * signed_pole);
            assert!((signed_frequency_entry(value, w, false) - expected).abs() < 1e-15);
            assert!((signed_frequency_entry(slope, w, true) - expected_slope).abs() < 1e-15);
        }
    }

    #[test]
    fn real_dc_and_imaginary_derivative_have_compatible_two_sided_limits() {
        let value = C64::from_re(0.9);
        let derivative = C64::new(0.0, -0.4);
        for w in [-f64::MIN_POSITIVE, -0.0, 0.0, f64::MIN_POSITIVE] {
            assert_eq!(signed_frequency_entry(value, w, false), value);
            assert_eq!(signed_frequency_entry(derivative, w, true), derivative);
        }
    }
}
