//! Element contractions, using the same P1 mass rules as conduction assembly.
//! These are control derivatives, not a replacement physical solver.

pub(super) fn source_mass(volume: f64, adjoint: [f64; 4]) -> [f64; 4] {
    let sum = adjoint.iter().sum::<f64>();
    adjoint.map(|value| (volume / 20.0) * (sum + value))
}

/// (dJ/dh, dJ/dT_ref, dJ/dq_out) on one uniform triangular face.
pub(super) fn boundary(
    area: f64, adjoint: [f64; 3], temperature: [f64; 3], h: f64, reference: f64,
) -> [f64; 3] {
    let sum = adjoint.iter().sum::<f64>();
    let delta = temperature.map(|value| value - reference);
    let diagonal = adjoint.iter().zip(delta).map(|(a, d)| a * d).sum::<f64>();
    let h_bar = -(area / 12.0) * (sum * delta.iter().sum::<f64>() + diagonal);
    [h_bar, h * (area / 3.0) * sum, -(area / 3.0) * sum]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn source_transpose_matches_independently_expanded_mass() {
        let lambda = [0.4, -2.0, 3.0, 0.125];
        let source = [2.0, 7.0, -1.0, 4.0];
        let volume = 3.5;
        let mut forward = 0.0;
        for i in 0..4 {
            for j in 0..4 {
                forward += lambda[i] * source[j] * volume / 20.0
                    * if i == j { 2.0 } else { 1.0 };
            }
        }
        let reverse = source_mass(volume, lambda).iter().zip(source)
            .map(|(a, f)| a * f).sum::<f64>();
        assert!((forward - reverse).abs() < 1e-12);
    }

    #[test]
    fn uniform_generation_uses_quarter_volume_not_a_diagonal_mass() {
        let values = source_mass(8.0, [1.0, 2.0, 3.0, 4.0]);
        assert!((values.iter().sum::<f64>() - 20.0).abs() < 1e-14);
    }

    #[test]
    fn robin_derivatives_include_consistent_off_diagonal_face_mass() {
        let a = [1.0, -2.0, 3.0];
        let t = [310.0, 320.0, 330.0];
        let [h, reference, outward] = boundary(6.0, a, t, 4.0, 300.0);
        let mut expected = 0.0;
        for i in 0..3 {
            for j in 0..3 {
                expected -= a[i] * (t[j] - 300.0) * 6.0 / 12.0
                    * if i == j { 2.0 } else { 1.0 };
            }
        }
        assert_eq!(h, expected);
        assert_eq!(reference, 16.0);
        assert_eq!(outward, -4.0);
    }

    #[test]
    fn equilibrium_has_zero_htc_derivative_but_not_zero_ambient_derivative() {
        let [h, reference, outward] = boundary(3.0, [1.0; 3], [300.0; 3], 7.0, 300.0);
        assert_eq!(h, 0.0);
        assert_eq!(reference, 21.0);
        assert_eq!(outward, -3.0);
    }
}
