//! Bounded linear face traces on the actual nonuniform axial geometry.
//!
//! Both choices use the same half-cell distance at the two faces. Limiting
//! that half jump by each neighboring centre jump keeps both traces in their
//! local intervals. The existing half-Courant transport bound therefore also
//! applies to the less dissipative centered choice.

/// Spatial limiter for native axial enthalpy transport.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum AxialLimiter {
    /// Smaller one-sided physical slope; the original default arithmetic.
    #[default]
    Minmod,
    /// Nonuniform centered slope limited by both neighboring centre values.
    /// On equal cells this is the monotonized-central (MC) limiter. Smooth
    /// fronts retain more structure than minmod without new face extrema.
    MonotonizedCentral,
}

/// Half-cell enthalpy increment and its derivatives in previous/centre/next
/// enthalpy coordinates. Inputs are finite jumps and geometric fractions
/// `half_width / centre_distance`, at most one. An extreme length ratio may
/// round a fraction to zero.
pub(super) fn half_slope(
    limiter: AxialLimiter,
    left: f64,
    right: f64,
    left_fraction: f64,
    right_fraction: f64,
) -> (f64, [f64; 3]) {
    if !((left > 0.0 && right > 0.0) || (left < 0.0 && right < 0.0)) {
        return (0.0, [0.0; 3]);
    }
    if limiter == AxialLimiter::Minmod {
        let a = left * left_fraction;
        let b = right * right_fraction;
        return if a.abs() <= b.abs() {
            (a, [-left_fraction, left_fraction, 0.0])
        } else {
            (b, [0.0, -right_fraction, right_fraction])
        };
    }

    // The derivative of the quadratic through the three actual centres is
    // (d_right*s_left + d_left*s_right)/(d_left+d_right). Work with bounded
    // fractions instead of dividing enthalpy by very short physical lengths.
    let sum = left_fraction + right_fraction;
    if sum == 0.0 {
        // Both physical gradients are below geometric representability.
        return (0.0, [0.0; 3]);
    }
    let a = (left_fraction / sum) * left_fraction;
    let b = (right_fraction / sum) * right_fraction;
    let centered = a * left + b * right;
    // On nonuniform grids, 2*minmod(s_left,s_right) alone is insufficient:
    // a large cell's face may be almost a whole centre distance away. These
    // centre-jump bounds retain monotone traces even at extreme grid ratios.
    if left.abs() <= centered.abs() && left.abs() <= right.abs() {
        (left, [-1.0, 1.0, 0.0])
    } else if right.abs() <= centered.abs() {
        (right, [0.0, -1.0, 1.0])
    } else {
        (centered, [-a, a - b, b])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nonuniform_affine_fields_and_both_face_bounds() {
        for lengths in [[1.0, 1.0, 1.0], [0.01, 10.0, 0.03], [10.0, 0.01, 20.0]] {
            let [previous, centre, next] = lengths;
            let dl = 0.5 * (previous + centre);
            let dr = 0.5 * (centre + next);
            for gradient in [-3.0, 0.0, 4.0] {
                let expected = gradient * 0.5 * centre;
                for limiter in [AxialLimiter::Minmod, AxialLimiter::MonotonizedCentral] {
                    let (slope, _) = half_slope(
                        limiter,
                        gradient * dl,
                        gradient * dr,
                        0.5 * centre / dl,
                        0.5 * centre / dr,
                    );
                    assert!((slope - expected).abs() <= 8e-15 * expected.abs().max(1.0));
                }
            }
            for left in [-100.0_f64, -0.01, 0.0, 0.01, 100.0] {
                for right in [-100.0_f64, -0.01, 0.0, 0.01, 100.0] {
                    let (slope, _) = half_slope(
                        AxialLimiter::MonotonizedCentral,
                        left,
                        right,
                        0.5 * centre / dl,
                        0.5 * centre / dr,
                    );
                    assert!(slope.abs() <= left.abs().min(right.abs()));
                    if left * right <= 0.0 {
                        assert_eq!(slope, 0.0);
                    }
                }
            }
        }
    }

    #[test]
    fn active_tangents_match_independent_perturbations() {
        for values in [
            [1.0, 2.0, 3.5],
            [1.0, 1.1, 9.0],
            [1.0, 8.9, 9.0],
            [9.0, 8.0, 6.5],
        ] {
            for fractions in [[0.5, 0.5], [0.8, 0.3]] {
                for limiter in [AxialLimiter::Minmod, AxialLimiter::MonotonizedCentral] {
                    let evaluate = |h: [f64; 3]| {
                        half_slope(
                            limiter,
                            h[1] - h[0],
                            h[2] - h[1],
                            fractions[0],
                            fractions[1],
                        )
                    };
                    let (_, derivative) = evaluate(values);
                    assert!(derivative.iter().sum::<f64>().abs() < 1e-15);
                    for axis in 0..3 {
                        let mut plus = values;
                        let mut minus = values;
                        plus[axis] += 1e-6;
                        minus[axis] -= 1e-6;
                        let numerical = (evaluate(plus).0 - evaluate(minus).0) / 2e-6;
                        assert!((numerical - derivative[axis]).abs() < 2e-9);
                    }
                }
            }
        }
    }

    #[test]
    fn large_same_sign_jumps_need_no_gradient_overflow() {
        for sign in [-1.0, 1.0] {
            let (slope, _) = half_slope(
                AxialLimiter::MonotonizedCentral,
                sign * f64::MAX,
                sign * f64::MAX,
                0.99,
                0.99,
            );
            assert!(slope.is_finite());
            assert_eq!(slope.signum(), sign);
        }
    }

    #[test]
    fn underflowed_geometric_fractions_retain_finite_donor_traces() {
        for fractions in [[0.0, 0.0], [0.0, 0.5], [0.5, 0.0]] {
            let (slope, tangent) = half_slope(
                AxialLimiter::MonotonizedCentral,
                1.0,
                2.0,
                fractions[0],
                fractions[1],
            );
            assert!(slope.is_finite());
            assert!(tangent.iter().all(|entry| entry.is_finite()));
            if fractions == [0.0, 0.0] {
                assert_eq!((slope, tangent), (0.0, [0.0; 3]));
            }
        }
    }

    #[test]
    fn periodic_heat_pulse_is_conservative_bounded_and_less_smeared() {
        // Independent finite-volume SSPRK2 evolution exercises the production
        // face kernel against an exactly translated cell-average pulse.
        let n = 64;
        let initial: Vec<_> = (0..n)
            .map(|i| if i < 16 { 100.0 } else { 10.0 })
            .collect();
        for sign in [-1.0, 1.0] {
            let mut errors = [0.0; 2];
            for (index, limiter) in [AxialLimiter::Minmod, AxialLimiter::MonotonizedCentral]
                .into_iter()
                .enumerate()
            {
                let euler = |h: &[f64]| {
                    let slopes: Vec<_> = (0..n)
                        .map(|i| {
                            half_slope(
                                limiter,
                                h[i] - h[(i + n - 1) % n],
                                h[(i + 1) % n] - h[i],
                                0.5,
                                0.5,
                            )
                            .0
                        })
                        .collect();
                    let flux: Vec<_> = (0..n)
                        .map(|i| {
                            if sign > 0.0 {
                                h[i] + slopes[i]
                            } else {
                                let next = (i + 1) % n;
                                -(h[next] - slopes[next])
                            }
                        })
                        .collect();
                    (0..n)
                        .map(|i| h[i] - 0.25 * (flux[i] - flux[(i + n - 1) % n]))
                        .collect::<Vec<_>>()
                };
                let mut h = initial.clone();
                for _ in 0..64 {
                    let second = euler(&euler(&h));
                    for i in 0..n {
                        h[i] = 0.5 * h[i] + 0.5 * second[i];
                        assert!((10.0 - 1e-12..=100.0 + 1e-12).contains(&h[i]));
                    }
                    assert!((h.iter().sum::<f64>() - initial.iter().sum::<f64>()).abs() < 1e-9);
                }
                errors[index] = (0..n)
                    .map(|i| {
                        let origin = if sign > 0.0 {
                            (i + n - 16) % n
                        } else {
                            (i + 16) % n
                        };
                        (h[i] - initial[origin]).abs()
                    })
                    .sum::<f64>()
                    / n as f64;
            }
            assert!(
                errors[1] < 0.85 * errors[0],
                "minmod={}, MC={}",
                errors[0],
                errors[1]
            );
        }
    }
}
