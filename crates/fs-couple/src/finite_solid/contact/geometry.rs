//! Relative material witnesses before rounding their absolute world positions.

use fs_math::dd::Dd;

/// `(first_origin + first_offset) - (second_origin + second_offset)` in one
/// inertial frame. Origins are retained coordinates; offsets are rotated rigid
/// arms or native nodal displacements. Rounding each world point first can erase
/// a real contact gap or lateral offset under a common translation and thereby
/// activate an unphysical impulse. No geometry is moved or tolerance enlarged.
pub(crate) fn separation(
    first_origin: [f64; 3],
    first_offset: [f64; 3],
    second_origin: [f64; 3],
    second_offset: [f64; 3],
) -> Option<[f64; 3]> {
    let mut result = [0.0; 3];
    for i in 0..3 {
        let value = (Dd::from_f64(first_origin[i]) - Dd::from_f64(second_origin[i]))
            + (Dd::from_f64(first_offset[i]) - Dd::from_f64(second_offset[i]));
        if !value.hi.is_finite() || !value.lo.is_finite() {
            return None;
        }
        result[i] = value.to_f64();
        if !result[i].is_finite() {
            return None;
        }
    }
    Some(result)
}

#[cfg(test)]
mod tests {
    use super::separation;

    #[test]
    fn g3_common_translation_cannot_erase_contact_or_lateral_offsets() {
        let a = [0.375, -0.125, 0.25];
        let b = [0.25, 0.125, 0.25];
        let expected = [0.125, -0.25, 0.0];
        for shift in [0.0, 2.0_f64.powi(53), -2.0_f64.powi(53)] {
            let origin = [shift; 3];
            assert_eq!(separation(origin, a, origin, b), Some(expected));
            assert_eq!(separation(origin, b, origin, a), Some(expected.map(|v| -v)));
            if shift != 0.0 {
                let rounded: [f64; 3] = core::array::from_fn(|i| (shift + a[i]) - (shift + b[i]));
                assert_eq!(
                    rounded, [0.0; 3],
                    "the previous world-point subtraction loses the gap"
                );
            }
        }
    }

    #[test]
    fn g0_local_offsets_survive_cancellation_of_distinct_origins() {
        let origin = 2.0_f64.powi(53);
        let a = [origin + 2.0, origin - 4.0, origin];
        let b = [origin; 3];
        assert_eq!(
            separation(a, [-1.875, 3.75, -0.125], b, [0.0; 3]),
            Some([0.125, -0.25, -0.125]),
        );
        // A native solid displacement is a retained offset, just like a rigid arm.
        assert_eq!(
            separation(b, [0.5, 0.0, 0.0], b, [0.375, 0.0, 0.0]),
            Some([0.125, 0.0, 0.0]),
        );
    }

    #[test]
    fn g0_nonfinite_or_unrepresentable_separation_refuses() {
        for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(separation([value, 0.0, 0.0], [0.0; 3], [0.0; 3], [0.0; 3]).is_none());
            assert!(separation([0.0; 3], [value, 0.0, 0.0], [0.0; 3], [0.0; 3]).is_none());
        }
        assert!(separation([f64::MAX; 3], [0.0; 3], [-f64::MAX; 3], [0.0; 3]).is_none());
    }
}
