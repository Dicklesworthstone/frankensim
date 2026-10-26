//! Stable convexity loss for an admitted radial, piecewise-affine H(B) law.
//! The source owner supplies its original segment slopes. Every contribution
//! is nonnegative; subtracting two large stored energies is unnecessary.

pub(super) fn bregman_density(
    before: [f64; 2],
    after: [f64; 2],
    h_after: f64,
    segments: impl IntoIterator<Item = (f64, f64, f64)>,
) -> Option<f64> {
    let a = before[0].hypot(before[1]);
    let b = after[0].hypot(after[1]);
    if !a.is_finite() || !b.is_finite() || !h_after.is_finite() || h_after < 0.0 {
        return None;
    }
    let low = a.min(b);
    let high = a.max(b);
    let mut radial = 0.0;
    for (left, right, slope) in segments {
        let lo = low.max(left);
        let hi = high.min(right);
        if hi <= lo {
            continue;
        }
        // Integral |r-a| H'(r) dr over this source segment. Factored widths
        // avoid cancellation in r^2-a^2 for nearly identical strong fields.
        let weight = 0.5 * (lo - a).abs() + 0.5 * (hi - a).abs();
        let term = slope * (hi - lo) * weight;
        if !term.is_finite() || term < 0.0 {
            return None;
        }
        radial += term;
    }
    let angular = if a == 0.0 || b == 0.0 {
        0.0
    } else {
        // |u-v|^2/2 = 1-cos(theta), stable even for a tiny rotation.
        let dx = before[0] / a - after[0] / b;
        let dy = before[1] / a - after[1] / b;
        (h_after * a) * (0.5 * dx.mul_add(dx, dy * dy))
    };
    let result = radial + angular;
    (result.is_finite() && result >= 0.0).then_some(result)
}

#[cfg(test)]
mod tests {
    use super::bregman_density;

    fn loss(before: [f64; 2], after: [f64; 2]) -> f64 {
        let b = after[0].hypot(after[1]);
        let h = if b <= 1.0 { b } else { 1.0 + 4.0 * (b - 1.0) };
        bregman_density(before, after, h, [(0.0, 1.0, 1.0), (1.0, 4.0, 4.0)]).unwrap()
    }
    fn energy(b: f64) -> f64 {
        if b <= 1.0 {
            0.5 * b * b
        } else {
            0.5 + (b - 1.0) + 2.0 * (b - 1.0).powi(2)
        }
    }
    fn near(a: f64, b: f64) {
        assert!((a - b).abs() <= 2e-13 * (1.0 + b.abs()), "{a} != {b}");
    }

    #[test]
    fn independent_piecewise_energy_identity_in_both_directions_and_rotations() {
        for before in [[0.0, 0.0], [0.5, 0.0], [2.0, 0.0], [-1.2, 0.5]] {
            for after in [[0.0, 0.0], [0.5, 0.0], [2.0, 0.0], [0.0, -2.5]] {
                let a = f64::hypot(before[0], before[1]);
                let b = f64::hypot(after[0], after[1]);
                let h = if b <= 1.0 { b } else { 1.0 + 4.0 * (b - 1.0) };
                let work = if b == 0.0 {
                    0.0
                } else {
                    h * ((after[0] / b) * (after[0] - before[0])
                        + (after[1] / b) * (after[1] - before[1]))
                };
                near(loss(before, after), work - energy(b) + energy(a));
            }
        }
        near(loss([0.5, 0.0], [2.0, 0.0]), 4.125);
        near(loss([2.0, 0.0], [0.5, 0.0]), 2.625);
    }

    #[test]
    fn linear_limit_is_half_the_squared_vector_increment_over_permeability() {
        let before = [1.2, -0.7];
        let after = [-0.4, 1.8];
        let reluctivity = 3.0;
        let actual = bregman_density(
            before,
            after,
            reluctivity * f64::hypot(after[0], after[1]),
            [(0.0, 4.0, reluctivity)],
        )
        .unwrap();
        near(
            actual,
            0.5 * reluctivity * ((after[0] - before[0]).powi(2) + (after[1] - before[1]).powi(2)),
        );
        assert_eq!(loss(before, before), 0.0);
    }

    #[test]
    fn tiny_steps_in_strong_fields_and_small_rotations_keep_positive_loss() {
        let before = [1e8, 0.0];
        let after = [1e8 + 1.0, 0.0];
        let actual = bregman_density(before, after, 3.0 * after[0], [(0.0, 2e8, 3.0)]).unwrap();
        assert_eq!(actual, 1.5);
        let rotation = loss([1.0, 1e-8], [1.0, -1e-8]);
        assert!((rotation / 2e-16 - 1.0).abs() < 1e-14);
    }

    #[test]
    fn unrepresentable_loss_is_refused() {
        assert!(bregman_density([f64::INFINITY, 0.0], [0.0, 0.0], 0.0, []).is_none());
        assert!(bregman_density([1.0, 0.0], [-1.0, 0.0], f64::MAX, [(0.0, 2.0, 1.0)]).is_none());
    }
}
