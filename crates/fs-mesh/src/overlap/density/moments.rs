//! Positive cubic simplex integration for a linearly varying density.
//!
//! Coordinates are relative to the caller's moment origin. Geometry clipping
//! remains in the parent; this allocation-free leaf integrates one tetrahedron.

/// Volume, mass, first mass moment and row-major second mass moment.
/// Positive relative densities are multiplied by the original density scale.
pub(crate) fn integrate(
    volume: f64,
    positions: [[f64; 3]; 4],
    relative_density: [f64; 4],
    density_scale: f64,
) -> Option<[f64; 14]> {
    if !volume.is_finite()
        || volume <= 0.0
        || !density_scale.is_finite()
        || density_scale <= 0.0
        || relative_density
            .iter()
            .any(|r| !r.is_finite() || *r <= 0.0 || *r > 1.0)
        || positions.iter().flatten().any(|x| !x.is_finite())
    {
        return None;
    }
    let sum_density: f64 = relative_density.iter().sum();
    let mass = product([volume, density_scale, sum_density * 0.25, 1.0]);
    if !mass.is_finite() || mass <= 0.0 {
        return None;
    }
    let mut out = [0.0; 14];
    out[0] = volume;
    out[1] = mass;
    let scale = positions
        .iter()
        .flatten()
        .fold(0.0_f64, |a, x| a.max(x.abs()));
    if scale == 0.0 {
        return Some(out);
    }
    let r = positions.map(|p| p.map(|x| x / scale));
    let sum: [f64; 3] = core::array::from_fn(|j| r.iter().map(|p| p[j]).sum());
    let shares = relative_density.map(|rho| rho / sum_density);
    for j in 0..3 {
        let mean: f64 = (0..4).map(|a| shares[a] * ((sum[j] + r[a][j]) / 5.0)).sum();
        out[2 + j] = product([mass, scale, mean, 1.0]);
        for k in j..3 {
            let square_sum: f64 = r.iter().map(|p| p[j] * p[k]).sum();
            // Integral(N_a r_j r_k) = V/120 times this bracket.
            // Dividing by M = V sum(rho)/4 gives the factor 1/30.
            // For j == k every summand is nonnegative, even when the origin
            // lies inside the tetrahedron; no raw-moment subtraction occurs.
            let moment: f64 = (0..4)
                .map(|a| {
                    shares[a]
                        * (((sum[j] + r[a][j]) * (sum[k] + r[a][k])
                            + square_sum
                            + r[a][j] * r[a][k])
                            / 30.0)
                })
                .sum();
            let value = product([mass, scale, scale, moment]);
            out[5 + 3 * j + k] = value;
            out[5 + 3 * k + j] = value;
        }
    }
    out.iter().all(|x| x.is_finite()).then_some(out)
}

// Binary factoring keeps finite products when an intermediate density*volume
// or squared coordinate would overflow/underflow. Inputs are already finite.
// Only the final result rounds into the subnormal range.
fn product(factors: [f64; 4]) -> f64 {
    let mut sign = 1.0;
    let mut fraction = 1.0;
    let mut exponent = 0_i32;
    for value in factors {
        if value == 0.0 {
            return 0.0;
        }
        if value < 0.0 {
            sign = -sign;
        }
        let mut value = value.abs();
        let mut shift = 0;
        if value < f64::MIN_POSITIVE {
            value *= 18_014_398_509_481_984.0;
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
    let value = if exponent > 1023 {
        f64::INFINITY
    } else if exponent >= -1022 {
        fraction * f64::from_bits(((exponent + 1023) as u64) << 52)
    } else if exponent < -1075 {
        0.0
    } else {
        let scale = f64::from_bits(((exponent + 2045) as u64) << 52);
        (fraction * scale) * f64::MIN_POSITIVE
    };
    sign * value
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64) {
        assert!(
            (a - b).abs() <= 4e-14 * a.abs().max(b.abs()).max(1e-300),
            "{a:e} != {b:e}"
        );
    }

    fn monomial(volume: f64, indices: &[usize]) -> f64 {
        let mut powers = [0; 4];
        for &i in indices {
            powers[i] += 1;
        }
        let factorial = |n: usize| (1..=n).fold(1.0, |v, x| v * x as f64);
        6.0 * volume * powers.map(factorial).iter().product::<f64>() / factorial(3 + indices.len())
    }

    #[test]
    fn graded_mass_moments_match_independent_cubic_integrals() {
        let positions = [
            [0.3, -0.8, 1.1],
            [1.7, 0.2, -0.4],
            [-0.5, 2.1, 0.6],
            [0.8, -0.3, -1.4],
        ];
        let rho = [2.0, 5.0, 3.0, 4.0];
        let volume = 0.27;
        let actual = integrate(volume, positions, rho.map(|r| r / 5.0), 5.0).unwrap();
        close(actual[1], volume * 3.5);
        for j in 0..3 {
            let mut first = 0.0;
            for a in 0..4 {
                for b in 0..4 {
                    first += rho[a] * positions[b][j] * monomial(volume, &[a, b]);
                }
            }
            close(actual[2 + j], first);
            for k in 0..3 {
                let mut second = 0.0;
                for a in 0..4 {
                    for b in 0..4 {
                        for c in 0..4 {
                            second += rho[a]
                                * positions[b][j]
                                * positions[c][k]
                                * monomial(volume, &[a, b, c]);
                        }
                    }
                }
                close(actual[5 + 3 * j + k], second);
            }
        }
    }

    #[test]
    fn unit_tet_retains_density_weighted_centroid_and_rotational_energy() {
        let points = [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let m = integrate(1.0 / 6.0, points, [0.25, 0.5, 0.75, 1.0], 4.0).unwrap();
        close(m[1], 5.0 / 12.0);
        close(m[2], 1.0 / 10.0);
        close(m[3], 13.0 / 120.0);
        close(m[4], 7.0 / 60.0);
        close(m[5], 7.0 / 180.0);
        close(m[6], 1.0 / 48.0);
        // u = [-y, x, 0], so K = (I_xx + I_yy)/2.
        close(0.5 * (m[5] + m[9]), 1.0 / 24.0);
        assert_ne!(m[2] / m[1], 0.25);
    }

    #[test]
    fn vertex_permutation_and_density_rescaling_preserve_physical_moments() {
        let p = [
            [-1.0, 0.5, 2.0],
            [1.0, 0.0, 2.0],
            [0.0, 2.0, -1.0],
            [0.0, 0.0, 1.0],
        ];
        let r = [0.2, 0.4, 0.8, 1.0];
        let a = integrate(0.4, p, r, 10.0).unwrap();
        let order = [2, 0, 3, 1];
        let b = integrate(0.4, order.map(|i| p[i]), order.map(|i| r[i]), 10.0).unwrap();
        let c = integrate(0.4, p, r.map(|x| x * 0.5), 20.0).unwrap();
        for i in 0..14 {
            close(a[i], b[i]);
            close(a[i], c[i]);
        }
    }

    #[test]
    fn origin_translation_obeys_first_and_second_moment_transport() {
        let p = [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let shift = [2.0, -3.0, 0.7];
        let q = p.map(|point| core::array::from_fn(|j| point[j] - shift[j]));
        let a = integrate(0.2, p, [0.2, 0.4, 0.8, 1.0], 5.0).unwrap();
        let b = integrate(0.2, q, [0.2, 0.4, 0.8, 1.0], 5.0).unwrap();
        for j in 0..3 {
            close(b[2 + j], a[2 + j] - a[1] * shift[j]);
            for k in 0..3 {
                close(
                    b[5 + 3 * j + k],
                    a[5 + 3 * j + k] - shift[j] * a[2 + k] - shift[k] * a[2 + j]
                        + a[1] * shift[j] * shift[k],
                );
            }
        }
    }

    #[test]
    fn binary_factoring_preserves_representable_density_length_scales() {
        let high = integrate(1e-300, [[1e150, -1e150, 0.0]; 4], [1.0; 4], 1.0).unwrap();
        let low = integrate(1e300, [[1e-150, -1e-150, 0.0]; 4], [1.0; 4], 1.0).unwrap();
        for m in [high, low] {
            close(m[5], 1.0);
            close(m[6], -1.0);
            close(m[9], 1.0);
        }
        close(
            integrate(1e-200, [[1.0; 3]; 4], [1.0; 4], 1e200).unwrap()[1],
            1.0,
        );
        assert_eq!(product([f64::from_bits(1), 0.75, 1.0, 1.0]).to_bits(), 1);
        assert_eq!(product([f64::from_bits(1), 0.5, 1.0, 1.0]), 0.0);
    }

    #[test]
    fn invalid_density_and_unrepresentable_storage_refuse() {
        assert!(integrate(0.0, [[0.0; 3]; 4], [1.0; 4], 1.0).is_none());
        assert!(integrate(1.0, [[0.0; 3]; 4], [0.0, 1.0, 1.0, 1.0], 1.0).is_none());
        assert!(integrate(1.0, [[0.0; 3]; 4], [1.0; 4], f64::NAN).is_none());
        assert!(integrate(1e300, [[1e300; 3]; 4], [1.0; 4], 1.0).is_none());
    }
}
