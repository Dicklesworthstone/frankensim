//! Exact P1 reference-density derivatives of a tetrahedron's kinetic form.
//!
//! This numerical leaf has no material lookup, solver, allocation or mutable
//! state. The caller supplies the original volume and a dimensionless weight.

/// `weight * integral(N_a * |u_h|^2 dV)` for each original vertex.
pub(crate) fn moments(volume: f64, weight: f64, displacement: [[f64; 3]; 4]) -> Option<[f64; 4]> {
    if !volume.is_finite()
        || volume <= 0.0
        || !weight.is_finite()
        || weight < 0.0
        || displacement.iter().flatten().any(|v| !v.is_finite())
    {
        return None;
    }
    let scale = displacement
        .iter()
        .flatten()
        .fold(0.0_f64, |a, v| a.max(v.abs()));
    if scale == 0.0 || weight == 0.0 {
        return Some([0.0; 4]);
    }
    let u = displacement.map(|v| v.map(|x| x / scale));
    let mut sum = [0.0; 3];
    let mut squares = 0.0;
    for v in u {
        for axis in 0..3 {
            sum[axis] += v[axis];
            squares = v[axis].mul_add(v[axis], squares);
        }
    }
    // Cubic simplex moments give the positive identity
    // integral(N_a |u_h|²) = V/120 * (|sum(u)+u_a|² + sum|u|² + |u_a|²).
    // Scaling before squaring retains finite results for large mode amplitudes
    // and small volumes/weights without constructing a 12x12 derivative tile.
    let result = u.map(|v| {
        let mut value = squares;
        for axis in 0..3 {
            let shifted = sum[axis] + v[axis];
            value = shifted.mul_add(shifted, v[axis].mul_add(v[axis], value));
        }
        positive_product([volume, weight, scale, scale, value / 120.0])
    });
    result.iter().all(|v| v.is_finite()).then_some(result)
}

// All factors are finite and nonnegative. Normalize each one before multiplying
// so a representable answer is not lost to an intermediate square or V*weight.
fn positive_product(factors: [f64; 5]) -> f64 {
    let mut fraction = 1.0;
    let mut exponent = 0_i32;
    for mut value in factors {
        if value == 0.0 {
            return 0.0;
        }
        let mut shift = 0;
        if value < f64::MIN_POSITIVE {
            value *= 18_014_398_509_481_984.0; // 2^54, exact.
            shift = -54;
        }
        let bits = value.to_bits();
        exponent += ((bits >> 52) & 0x7ff) as i32 - 1023 + shift;
        fraction *= f64::from_bits((bits & 0x000f_ffff_ffff_ffff) | (1023_u64 << 52));
        if fraction >= 2.0 {
            fraction *= 0.5;
            exponent += 1;
        }
    }
    if exponent > 1023 {
        f64::INFINITY
    } else if exponent >= -1022 {
        fraction * f64::from_bits(((exponent + 1023) as u64) << 52)
    } else if exponent < -1075 {
        0.0
    } else {
        // One final rounding into the subnormal range, including the half-ULP
        // tie at exponent -1075. The first product is always normal.
        let scale = f64::from_bits(((exponent + 2045) as u64) << 52);
        (fraction * scale) * f64::MIN_POSITIVE
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64) {
        assert!(
            (a - b).abs() <= 3e-14 * a.abs().max(b.abs()).max(1e-300),
            "{a:e} != {b:e}"
        );
    }

    #[test]
    fn independent_cubic_simplex_integrals_match_all_four_density_derivatives() {
        let u = [
            [0.3, -0.8, 1.1],
            [1.7, 0.2, -0.4],
            [-0.5, 2.1, 0.6],
            [0.8, -0.3, -1.4],
        ];
        let volume = 0.27;
        let weight = 0.63;
        let actual = moments(volume, weight, u).unwrap();
        for a in 0..4 {
            let mut expected = 0.0;
            for i in 0..4 {
                for j in 0..4 {
                    // Independent Dirichlet integral: 6 V * product(alpha!) / 6!.
                    let mut powers = [0; 4];
                    powers[a] += 1;
                    powers[i] += 1;
                    powers[j] += 1;
                    let factorial = |n| match n {
                        0 | 1 => 1.0,
                        2 => 2.0,
                        3 => 6.0,
                        _ => unreachable!(),
                    };
                    let integral = volume / 120.0 * powers.map(factorial).iter().product::<f64>();
                    for axis in 0..3 {
                        expected += weight * integral * u[i][axis] * u[j][axis];
                    }
                }
            }
            close(actual[a], expected);
        }
    }

    #[test]
    fn constant_translation_and_cartesian_rotation_retain_physical_moments() {
        for value in moments(1.0 / 6.0, 0.7, [[2.0, -1.0, 3.0]; 4]).unwrap() {
            close(value, 0.7 * 14.0 / 24.0);
        }
        let rotation = [[0.0; 3], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0; 3]];
        let gradient = moments(1.0 / 6.0, 1.0, rotation).unwrap();
        let density = [2.0, 5.0, 3.0, 4.0];
        let inertia: f64 = gradient.into_iter().zip(density).map(|(g, r)| g * r).sum();
        close(inertia, 11.0 / 90.0);
        close(gradient[0], 1.0 / 180.0);
        close(gradient[1], 1.0 / 90.0);
    }

    #[test]
    fn original_vertex_permutation_preserves_the_corresponding_derivatives() {
        let u = [
            [1.0, 0.0, 2.0],
            [0.5, -2.0, 1.0],
            [-1.0, 0.3, 0.0],
            [0.0, 1.0, 1.5],
        ];
        let permutation = [2, 0, 3, 1];
        let a = moments(0.4, 0.8, u).unwrap();
        let b = moments(0.4, 0.8, permutation.map(|i| u[i])).unwrap();
        for i in 0..4 {
            close(b[i], a[permutation[i]]);
        }
    }

    #[test]
    fn factored_evaluation_retains_representable_extreme_density_modes() {
        let high = moments(1e-200, 1e-200, [[1e200, 0.0, 0.0]; 4]).unwrap();
        let low = moments(1e200, 1e200, [[1e-200, 0.0, 0.0]; 4]).unwrap();
        for value in high.into_iter().chain(low) {
            close(value, 0.25);
        }
        assert_eq!(
            positive_product([f64::from_bits(1), 1.0, 1.0, 1.0, 1.0]).to_bits(),
            1
        );
        assert_eq!(
            positive_product([f64::from_bits(1), 0.5, 1.0, 1.0, 1.0]),
            0.0
        );
        assert_eq!(
            positive_product([f64::from_bits(1), 0.75, 1.0, 1.0, 1.0]).to_bits(),
            1
        );
        assert!(moments(1e200, 1.0, [[1e200, 0.0, 0.0]; 4]).is_none());
    }

    #[test]
    fn zero_motion_zero_weight_and_nonfinite_inputs_have_explicit_results() {
        assert_eq!(moments(1.0, 0.0, [[3.0; 3]; 4]), Some([0.0; 4]));
        assert_eq!(moments(1.0, 1.0, [[0.0; 3]; 4]), Some([0.0; 4]));
        assert!(moments(0.0, 1.0, [[1.0; 3]; 4]).is_none());
        assert!(moments(1.0, -1.0, [[1.0; 3]; 4]).is_none());
        assert!(moments(1.0, 1.0, [[f64::NAN; 3]; 4]).is_none());
    }
}
