//! SI frame and finite-band admission shared by both gas inlet consumers.

pub(super) fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

pub(super) fn rotation_is_proper(r: [[f64; 3]; 3]) -> bool {
    if !r.iter().flatten().all(|x| x.is_finite()) {
        return false;
    }
    for i in 0..3 {
        for j in 0..3 {
            let target = if i == j { 1.0 } else { 0.0 };
            if (dot(r[i], r[j]) - target).abs() > 256.0 * f64::EPSILON {
                return false;
            }
        }
    }
    let determinant = r[0][0] * (r[1][1] * r[2][2] - r[1][2] * r[2][1])
        - r[0][1] * (r[1][0] * r[2][2] - r[1][2] * r[2][0])
        + r[0][2] * (r[1][0] * r[2][1] - r[1][1] * r[2][0]);
    (determinant - 1.0).abs() <= 512.0 * f64::EPSILON
}

pub(super) fn to_field(
    r: [[f64; 3]; 3],
    origin: [f64; 3],
    point: [f64; 3],
    displacement: [f64; 3],
) -> [f64; 3] {
    core::array::from_fn(|a| {
        origin[a]
            + (0..3)
                .map(|i| r[i][a] * (point[i] + displacement[i]))
                .sum::<f64>()
    })
}

pub(super) fn to_grid(r: [[f64; 3]; 3], velocity: [f64; 3]) -> [f64; 3] {
    r.map(|row| dot(row, velocity))
}

pub(super) fn half_edges(
    r: [[f64; 3]; 3],
    spacing: [f64; 3],
    normal_axis: usize,
    point: bool,
) -> [[f64; 3]; 3] {
    core::array::from_fn(|axis| {
        if point || axis == normal_axis {
            [0.0; 3]
        } else {
            r[axis].map(|x| x * (0.5 * spacing[axis]))
        }
    })
}

pub(super) fn speed(v: [f64; 3]) -> f64 {
    let scale = v.iter().fold(0.0_f64, |m, x| m.max(x.abs()));
    if scale == 0.0 {
        0.0
    } else {
        let unit = v.map(|x| x / scale);
        scale * fs_math::det::sqrt(dot(unit, unit))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    const R: [[f64; 3]; 3] = [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];

    #[test]
    fn translated_rotated_frozen_wave_preserves_phase_and_velocity() {
        // A field-frame wave k=(2,3,-1) observed at a moving Cartesian face.
        // Grid x corresponds to field -y, so convecting +grid x moves +field y.
        let p = to_field(R, [10.0, 20.0, 30.0], [4.0, 5.0, 6.0], [-2.0, 0.0, 0.0]);
        assert_eq!(p, [15.0, 18.0, 36.0]);
        assert_eq!(dot([2.0, 3.0, -1.0], p), 48.0);
        assert_eq!(to_grid(R, [2.0, 3.0, -1.0]), [-3.0, 2.0, -1.0]);
        assert!(rotation_is_proper(R));
        assert!(!rotation_is_proper([
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, -1.0]
        ]));
        let mut scaled = R;
        scaled[0][1] = -1.01;
        assert!(!rotation_is_proper(scaled));
    }

    #[test]
    fn rectangular_face_filter_integrates_physical_tangent_edges() {
        // x-normal 2-by-6 m face: after a quarter-turn its field edges are
        // +/-1 in x and +/-3 in z, with no averaging through the boundary.
        let edges = half_edges(R, [0.1, 2.0, 6.0], 0, false);
        assert_eq!(edges, [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 0.0, 3.0]]);
        let k = [0.7, -0.2, 0.3];
        let analytic = (0.7_f64.sin() / 0.7) * (0.9_f64.sin() / 0.9);
        let n = 160;
        let mut integral = 0.0;
        for i in 0..n {
            for j in 0..n {
                let x = -1.0 + (i as f64 + 0.5) * 2.0 / n as f64;
                let z = -3.0 + (j as f64 + 0.5) * 6.0 / n as f64;
                integral += dot(k, [x, 0.0, z]).cos() / (n * n) as f64;
            }
        }
        assert!((integral - analytic).abs() < 8e-6);
        assert_eq!(half_edges(R, [0.1, 2.0, 6.0], 0, true), [[0.0; 3]; 3]);
    }

    #[test]
    fn advective_speed_avoids_square_overflow_and_is_rotation_invariant() {
        for v in [
            [3.0, 4.0, 12.0],
            [3e200, 4e200, 12e200],
            [3e-200, 4e-200, 12e-200],
        ] {
            let expected = 13.0 * (v[0] / 3.0);
            assert!((speed(v) / expected - 1.0).abs() < 3e-15);
            assert_eq!(speed(v), speed(to_grid(R, v)));
        }
        assert_eq!(speed([0.0; 3]), 0.0);
    }
}
