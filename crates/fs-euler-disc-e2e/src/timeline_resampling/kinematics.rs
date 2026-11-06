//! Continuous speed bounds for the original Hermite/normalized-SLERP pose.
//!
//! Momenta are a different reconstructed channel. Their interpolation is not
//! the derivative of the pose and cannot bound acoustic arrival compression.
//! Bounds use the existing interval arithmetic owner; no trajectory is changed.

use fs_ivl::Interval;

pub(super) const LINEAR_ROTATION_DOT: f64 = 1.0 - 1.0e-12;

/// Share the exact double-cover choice and branch predicate with the sampler.
pub(super) fn shortest_rotation(
    first: [f64; 4],
    mut second: [f64; 4],
) -> ([f64; 4], f64) {
    let mut dot = first
        .iter()
        .zip(second)
        .map(|(first, second)| first * second)
        .sum::<f64>();
    if dot < 0.0 {
        for component in &mut second {
            *component = -*component;
        }
        dot = -dot;
    }
    (second, dot.clamp(-1.0, 1.0))
}

/// Fixed work and storage; alpha bounds enclose the selected span endpoints.
pub(super) fn pose_speed_bounds(
    position: [[f64; 3]; 2],
    velocity: [[f64; 3]; 2],
    orientation: [[f64; 4]; 2],
    duration_s: f64,
    alpha: [Interval; 2],
) -> Option<[f64; 2]> {
    if !(duration_s > 0.0 && duration_s.is_finite())
        || position
            .iter()
            .flatten()
            .chain(velocity.iter().flatten())
            .chain(orientation.iter().flatten())
            .any(|value| !value.is_finite())
        || alpha.iter().any(|value| value.lo() < 0.0 || value.hi() > 1.0)
        || alpha[0].lo() > alpha[1].hi()
    {
        return None;
    }
    let translation = translation_speed_bound(position, velocity, duration_s, alpha);
    let rotation = rotation_speed_bound(orientation, duration_s)?;
    [translation, rotation]
        .iter()
        .all(|value| value.is_finite() && *value >= 0.0)
        .then_some([translation, rotation])
}

fn translation_speed_bound(
    position: [[f64; 3]; 2],
    velocity: [[f64; 3]; 2],
    duration_s: f64,
    alpha: [Interval; 2],
) -> f64 {
    if position[0] == position[1] && velocity == [[0.0; 3]; 2] {
        return 0.0;
    }
    // The derivative of cubic Hermite position is this quadratic Bernstein
    // polynomial. In particular, its middle control need not be bounded by
    // either accepted endpoint velocity.
    let controls = [
        velocity[0].map(Interval::point),
        core::array::from_fn(|axis| {
            Interval::point(3.0)
                * (Interval::point(position[1][axis]) - Interval::point(position[0][axis]))
                / Interval::point(duration_s)
                - Interval::point(velocity[0][axis])
                - Interval::point(velocity[1][axis])
        }),
        velocity[1].map(Interval::point),
    ];
    let clipped = [
        blossom(controls, alpha[0], alpha[0]),
        blossom(controls, alpha[0], alpha[1]),
        blossom(controls, alpha[1], alpha[1]),
    ];
    // Four fixed subspans tighten the convex-hull bound cheaply. This is a
    // continuous enclosure, not four velocity samples or an adaptive solver.
    let mut maximum = 0.0_f64;
    for half in split(clipped) {
        for quarter in split(half) {
            for control in quarter {
                maximum = maximum.max(norm_upper(control));
            }
        }
    }
    maximum
}

fn blossom(
    controls: [[Interval; 3]; 3],
    first: Interval,
    second: Interval,
) -> [Interval; 3] {
    let one = Interval::point(1.0);
    core::array::from_fn(|axis| {
        let left = (one - first) * controls[0][axis] + first * controls[1][axis];
        let right = (one - first) * controls[1][axis] + first * controls[2][axis];
        (one - second) * left + second * right
    })
}

fn split(controls: [[Interval; 3]; 3]) -> [[[Interval; 3]; 3]; 2] {
    let half = Interval::point(0.5);
    let left = core::array::from_fn(|axis| (controls[0][axis] + controls[1][axis]) * half);
    let right = core::array::from_fn(|axis| (controls[1][axis] + controls[2][axis]) * half);
    let middle = core::array::from_fn(|axis| (left[axis] + right[axis]) * half);
    [[controls[0], left, middle], [middle, right, controls[2]]]
}

fn rotation_speed_bound(orientation: [[f64; 4]; 2], duration_s: f64) -> Option<f64> {
    let [first, second] = orientation;
    let (second, dot) = shortest_rotation(first, second);
    if first == second {
        return Some(0.0);
    }
    let first = first.map(Interval::point);
    let second = second.map(Interval::point);
    let delta: [Interval; 4] = core::array::from_fn(|i| second[i] - first[i]);
    // For a normalized raw quaternion q, angular speed is
    // 2*|q wedge dq/dt|/|q|^2. The wedge of q0 and (q1-q0) avoids
    // subtracting nearly equal unit-quaternion Gram determinants.
    let mut wedge = [Interval::point(0.0); 6];
    let mut entry = 0;
    for i in 0..4 {
        for j in i + 1..4 {
            wedge[entry] = first[i] * delta[j] - first[j] * delta[i];
            entry += 1;
        }
    }
    let wedge = Interval::point(norm_upper(wedge));
    let first_norm = norm_squared(first);
    let second_norm = norm_squared(second);
    let (minimum_norm, factor) = if dot > LINEAR_ROTATION_DOT {
        // q=(1-u)q0+u*q1. Its squared norm is the linear interpolation
        // of endpoint squared norms minus u(1-u)*|q1-q0|^2.
        let lower = Interval::point(first_norm.lo().min(second_norm.lo()))
            - Interval::point(0.25) * norm_squared(delta);
        (lower.lo(), Interval::point(1.0))
    } else {
        // Use the SAME represented angle and denominator as SLERP. Account
        // for the final quaternion normalization instead of assuming the
        // rounded endpoint quaternions are exactly unit/orthogonal.
        let angle = dot.acos();
        let denominator = Interval::point(angle.sin());
        if denominator.lo() <= 0.0 {
            return None;
        }
        let angle = Interval::point(angle);
        let sine = angle.sin();
        let cosine = angle.cos();
        let weight = Interval::new(0.0, angle.hi()).sin() / denominator;
        let weights_squared = Interval::new(0.0, weight.hi()) * Interval::new(0.0, weight.hi());
        let product = first
            .iter()
            .zip(second)
            .fold(Interval::point(0.0), |sum, (a, b)| sum + *a * b);
        let scale = sine / denominator;
        // sin((1-u)a)^2 + sin(u*a)^2
        //   + 2*cos(a)*sin((1-u)a)*sin(u*a) = sin(a)^2.
        let lower = scale * scale
            + weights_squared * (first_norm - Interval::point(1.0))
            + weights_squared * (second_norm - Interval::point(1.0))
            + Interval::point(2.0) * weights_squared * (product - cosine);
        (lower.lo(), angle * sine / (denominator * denominator))
    };
    if !(minimum_norm > 0.0 && minimum_norm.is_finite()) {
        return None;
    }
    Some(
        (Interval::point(2.0) * factor * wedge
            / (Interval::point(duration_s) * Interval::point(minimum_norm)))
        .hi(),
    )
}

fn norm_squared<const N: usize>(vector: [Interval; N]) -> Interval {
    vector.iter().fold(Interval::point(0.0), |sum, value| {
        let magnitude = value.abs();
        sum + magnitude * magnitude
    })
}

fn norm_upper<const N: usize>(vector: [Interval; N]) -> f64 {
    let scale = vector.iter().fold(0.0_f64, |a, b| a.max(b.abs().hi()));
    if scale == 0.0 || !scale.is_finite() {
        return scale;
    }
    let scale = Interval::point(scale);
    let scaled = vector.map(|value| value / scale);
    (norm_squared(scaled).sqrt() * scale).hi()
}

#[cfg(test)]
mod tests {
    use super::*;

    const IDENTITY: [f64; 4] = [1.0, 0.0, 0.0, 0.0];

    fn span(first: f64, last: f64) -> [Interval; 2] {
        [Interval::point(first), Interval::point(last)]
    }

    #[test]
    fn hermite_interior_speed_is_not_endpoint_momentum() {
        let full = pose_speed_bounds(
            [[0.0; 3], [2.0, 0.0, 0.0]],
            [[0.0; 3]; 2],
            [IDENTITY; 2],
            1.0,
            span(0.0, 1.0),
        )
        .unwrap();
        assert!(full[0] >= 3.0 && full[0] < 3.0 + 1.0e-12);
        assert_eq!(full[1], 0.0);
        for i in 0..=128 {
            let u = f64::from(i) / 128.0;
            assert!(12.0 * u * (1.0 - u) <= full[0]);
        }
        let clipped = pose_speed_bounds(
            [[0.0; 3], [2.0, 0.0, 0.0]],
            [[0.0; 3]; 2],
            [IDENTITY; 2],
            1.0,
            span(0.0, 0.125),
        )
        .unwrap();
        assert!(clipped[0] >= 1.3125 && clipped[0] < 1.3125 + 1.0e-12);
    }

    #[test]
    fn rotation_bound_tracks_slerp_nlerp_and_double_cover() {
        for angle in [1.0e-13_f64, 1.0e-6, 1.0e-4, 0.7, core::f64::consts::PI - 1.0e-10] {
            let end = [(0.5 * angle).cos(), (0.5 * angle).sin(), 0.0, 0.0];
            let mut expected = None;
            for sign in [1.0, -1.0] {
                let bound = pose_speed_bounds(
                    [[0.0; 3]; 2],
                    [[0.0; 3]; 2],
                    [IDENTITY, end.map(|value| sign * value)],
                    0.25,
                    span(0.0, 1.0),
                )
                .unwrap();
                assert_eq!(bound[0], 0.0);
                assert!(bound[1] >= 4.0 * angle * (1.0 - 1.0e-14), "{angle}: {bound:?}");
                assert!(bound[1] < 4.0 * angle * (1.0 + 1.0e-5), "{angle}: {bound:?}");
                if let Some(expected) = expected {
                    assert_eq!(bound, expected);
                }
                expected = Some(bound);
            }
        }
    }

    #[test]
    fn vector_hermite_bound_encloses_the_continuous_derivative() {
        let position = [[1.0, -2.0, 3.0], [4.0, 2.0, -1.0]];
        let velocity = [[10.0, -4.0, 5.0], [-7.0, 8.0, -3.0]];
        let bound = pose_speed_bounds(
            position, velocity, [IDENTITY; 2], 0.5, span(0.05, 0.8),
        )
        .unwrap()[0];
        let mut sampled = 0.0_f64;
        for i in 0..=512 {
            let u = 0.05 + 0.75 * f64::from(i) / 512.0;
            let v: [f64; 3] = core::array::from_fn(|axis| {
                (6.0 * u * u - 6.0 * u) * (position[0][axis] - position[1][axis]) / 0.5
                    + (3.0 * u * u - 4.0 * u + 1.0) * velocity[0][axis]
                    + (3.0 * u * u - 2.0 * u) * velocity[1][axis]
            });
            let speed = v.iter().map(|value| value * value).sum::<f64>().sqrt();
            sampled = sampled.max(speed);
            assert!(speed <= bound, "{speed} > {bound}");
        }
        assert!(bound < 1.25 * sampled);
    }

    #[test]
    fn exact_rest_and_unrepresentable_motion_have_explicit_results() {
        let q = [0.6, 0.0, 0.8, 0.0];
        assert_eq!(
            pose_speed_bounds(
                [[1.0, -2.0, 3.0]; 2], [[0.0; 3]; 2],
                [q, q.map(|value| -value)], 1.0, span(0.0, 1.0),
            ),
            Some([0.0, 0.0]),
        );
        assert!(pose_speed_bounds(
            [[-f64::MAX, 0.0, 0.0], [f64::MAX, 0.0, 0.0]],
            [[0.0; 3]; 2], [IDENTITY; 2], 1.0, span(0.0, 1.0),
        ).is_none());
        assert!(pose_speed_bounds(
            [[0.0; 3]; 2], [[0.0; 3]; 2], [IDENTITY; 2],
            0.0, span(0.0, 1.0),
        ).is_none());
    }
}
