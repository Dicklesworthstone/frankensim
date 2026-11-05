//! Cartesian derivatives of the outgoing Helmholtz Green function.
//!
//! For exp(-i omega t), G = exp(ikr)/(4 pi r). A symmetric integrated
//! Lighthill stress Q gives p = Q_ij d_i d_j G and
//! u_l = Q_ij d_l d_i d_j G / (i rho omega). Both trace and shear survive.
//! These are exact point-tensor transfers, including their near fields.

use fs_math::{c64::C64, det};

/// Physical tensor entries, without engineering-shear factors.
pub(super) const AXES: [(usize, usize); 6] = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)];

/// Rows p, ux, uy, uz; columns Qxx, Qyy, Qzz, Qxy, Qxz, Qyz [N m].
pub(super) fn transfer(
    rho: f64,
    sound_speed: f64,
    frequency: f64,
    displacement: [f64; 3],
) -> Option<[[C64; 6]; 4]> {
    let omega = core::f64::consts::TAU * frequency;
    let k = omega / sound_speed;
    let rho_omega = rho * omega;
    let r = det::hypot(
        det::hypot(displacement[0], displacement[1]),
        displacement[2],
    );
    let inv_r = 1.0 / r;
    if [rho, sound_speed, frequency, omega, k, rho_omega, r, inv_r]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
        || displacement.iter().any(|x| !x.is_finite())
        || !(k * r).is_finite()
    {
        return None;
    }
    let n = displacement.map(|x| x / r);
    let (sin, cos) = det::sin_cos(k * r);
    let green = C64::new(cos, sin).scale((0.25 / core::f64::consts::PI) / r);
    let inv_r2 = inv_r * inv_r;
    let inv_r3 = inv_r2 * inv_r;
    let k2 = k * k;
    // Hessian: G [a n_i n_j + b delta_ij].
    let a = C64::new(3.0 * inv_r2 - k2, -3.0 * k * inv_r);
    let b = C64::new(-inv_r2, k * inv_r);
    // Third derivative: G [c n_i n_j n_l + d (delta_ij n_l +
    // delta_il n_j + delta_jl n_i)]. Differentiating the directions is
    // essential: pressure-null directions can carry transverse velocity.
    let c = C64::new(6.0 * k2 * inv_r - 15.0 * inv_r3, 15.0 * k * inv_r2 - k2 * k);
    let d = a.scale(inv_r);
    let mut field = [[C64::ZERO; 6]; 4];
    for (entry, (i, j)) in AXES.into_iter().enumerate() {
        let diagonal = if i == j { 1.0 } else { 0.0 };
        let multiplicity = if i == j { 1.0 } else { 2.0 };
        field[0][entry] = (green * (a.scale(n[i] * n[j]) + b.scale(diagonal))).scale(multiplicity);
        for l in 0..3 {
            let angular =
                diagonal * n[l] + if i == l { n[j] } else { 0.0 } + if j == l { n[i] } else { 0.0 };
            let gradient =
                (green * (c.scale(n[i] * n[j] * n[l]) + d.scale(angular))).scale(multiplicity);
            field[l + 1][entry] = C64::new(gradient.im / rho_omega, -gradient.re / rho_omega);
        }
    }
    field
        .iter()
        .flatten()
        .all(|x| x.re.is_finite() && x.im.is_finite())
        .then_some(field)
}

/// Instantaneous resolved rho u_i u_j dV, before any temporal transform.
pub(super) fn reynolds_moment(rho: f64, volume: f64, velocity: [f64; 3]) -> Option<[f64; 6]> {
    let mass = rho * volume;
    if [rho, volume, mass]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
        || velocity.iter().any(|x| !x.is_finite())
    {
        return None;
    }
    let moment = AXES.map(|(i, j)| (mass * velocity[i]) * velocity[j]);
    moment.iter().all(|x| x.is_finite()).then_some(moment)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn apply(x: [f64; 3], q: [C64; 6]) -> [C64; 4] {
        transfer(1.2, 343.0, 127.0, x).unwrap().map(|row| {
            row.into_iter()
                .zip(q)
                .fold(C64::ZERO, |s, (h, q)| s + h * q)
        })
    }

    fn close(a: C64, b: C64, relative: f64) {
        assert!(
            (a - b).abs() <= relative * (1e-12 + b.abs()),
            "{a:?} != {b:?}"
        );
    }

    fn green(x: [f64; 3]) -> C64 {
        let r = det::sqrt(x.iter().map(|x| x * x).sum());
        let phase = core::f64::consts::TAU * 127.0 * r / 343.0;
        let (s, c) = det::sin_cos(phase);
        C64::new(c, s).scale(1.0 / (4.0 * core::f64::consts::PI * r))
    }

    fn stress() -> [C64; 6] {
        [
            C64::new(0.7, 0.1),
            C64::new(-0.3, 0.2),
            C64::new(0.9, -0.4),
            C64::new(0.23, -0.15),
            C64::new(-0.11, 0.3),
            C64::new(0.4, 0.2),
        ]
    }

    #[test]
    fn g1_pressure_is_the_full_green_hessian_under_spatial_refinement() {
        let x = [1.3, -0.8, 0.6];
        let q = stress();
        let exact = apply(x, q)[0];
        let differences = [0.002, 0.001, 0.0005].map(|h| {
            let mut value = C64::ZERO;
            for (entry, (i, j)) in AXES.into_iter().enumerate() {
                let derivative = if i == j {
                    let mut plus = x;
                    plus[i] += h;
                    let mut minus = x;
                    minus[i] -= h;
                    (green(plus) - green(x).scale(2.0) + green(minus)).scale(1.0 / (h * h))
                } else {
                    let mut mixed = C64::ZERO;
                    for a in [-1.0, 1.0] {
                        for b in [-1.0, 1.0] {
                            let mut p = x;
                            p[i] += a * h;
                            p[j] += b * h;
                            mixed = mixed + green(p).scale(a * b);
                        }
                    }
                    mixed.scale(0.5 / (h * h))
                };
                value = value + derivative * q[entry];
            }
            (value - exact).abs()
        });
        assert!(differences[1] < 0.27 * differences[0], "{differences:?}");
        assert!(differences[2] < 0.28 * differences[1], "{differences:?}");
        assert!(
            differences[2] < 2e-6 * exact.abs(),
            "{differences:?}, exact {exact:?}"
        );
    }

    #[test]
    fn g1_velocity_obeys_all_three_components_of_linear_momentum() {
        let x = [1.3, -0.8, 0.6];
        let q = stress();
        let field = apply(x, q);
        let h = 1e-5;
        for axis in 0..3 {
            let mut plus = x;
            plus[axis] += h;
            let mut minus = x;
            minus[axis] -= h;
            let gradient = (apply(plus, q)[0] - apply(minus, q)[0]).scale(0.5 / h);
            close(
                gradient,
                C64::new(0.0, 1.2 * core::f64::consts::TAU * 127.0) * field[axis + 1],
                2e-9,
            );
        }
    }

    #[test]
    fn g1_isotropic_trace_retains_the_helmholtz_source() {
        let x = [1.3, -0.8, 0.6];
        let q = C64::new(0.7, -0.3);
        let field = apply(x, [q, q, q, C64::ZERO, C64::ZERO, C64::ZERO]);
        let k = core::f64::consts::TAU * 127.0 / 343.0;
        close(field[0], (green(x) * q).scale(-k * k), 3e-14);
    }

    #[test]
    fn g1_pressure_null_still_has_transverse_near_field_velocity() {
        let mut q = [C64::ZERO; 6];
        q[3] = C64::ONE;
        let field = apply([0.7, 0.0, 0.0], q);
        assert_eq!(field[0], C64::ZERO);
        assert_eq!(field[1], C64::ZERO);
        assert!(field[2].abs() > 1e-4);
        assert_eq!(field[3], C64::ZERO);
    }

    #[test]
    fn g3_cartesian_rotation_preserves_pressure_and_rotates_velocity() {
        let x = [1.3, -0.8, 0.6];
        let q = stress();
        let field = apply(x, q);
        // Proper cyclic rotation: (x,y,z) -> (z,x,y), including off-diagonal stress.
        let rotated = apply([x[2], x[0], x[1]], [q[2], q[0], q[1], q[4], q[5], q[3]]);
        close(rotated[0], field[0], 3e-15);
        for (actual, expected) in rotated[1..].iter().zip([field[3], field[1], field[2]]) {
            close(*actual, expected, 3e-15);
        }
    }

    #[test]
    fn g1_far_field_has_quadrupole_directivity_and_outgoing_impedance() {
        let q = stress();
        let k = core::f64::consts::TAU * 127.0 / 343.0;
        let x = [6000.0, 8000.0, 0.0];
        let n = [0.6, 0.8, 0.0];
        let field = apply(x, q);
        let projection = q[0].scale(0.36) + q[1].scale(0.64) + q[3].scale(0.96);
        close(field[0], (green(x) * projection).scale(-k * k), 2e-4);
        let radial = field[1].scale(n[0]) + field[2].scale(n[1]);
        close(radial.scale(1.2 * 343.0), field[0], 8e-5);
    }

    #[test]
    fn g1_resolved_velocity_is_squared_before_forming_the_source_spectrum() {
        let rho = 1.2;
        let volume = 0.004;
        let amplitude = 3.0;
        let mut first = C64::ZERO;
        let mut second = C64::ZERO;
        let mut mean = 0.0;
        for i in 0..64 {
            let t = core::f64::consts::TAU * i as f64 / 64.0;
            let moment = reynolds_moment(rho, volume, [amplitude * det::cos(t), 0.0, 0.0]).unwrap();
            first = first + C64::new(det::cos(t), det::sin(t)).scale(moment[0] / 32.0);
            second =
                second + C64::new(det::cos(2.0 * t), det::sin(2.0 * t)).scale(moment[0] / 32.0);
            mean += moment[0] / 64.0;
        }
        assert!(first.abs() < 1e-16);
        let expected = 0.5 * rho * volume * amplitude * amplitude;
        close(second, C64::from_re(expected), 2e-14);
        assert!((mean - expected).abs() < 1e-16);
    }

    #[test]
    fn g0_source_units_and_invalid_domains() {
        assert_eq!(
            reynolds_moment(2.0, 0.5, [1.0, -2.0, 3.0]),
            Some([1.0, 4.0, 9.0, -2.0, 3.0, -6.0])
        );
        for x in [0.0, -1.0, f64::INFINITY, f64::NAN] {
            assert!(reynolds_moment(x, 0.5, [1.0; 3]).is_none());
            assert!(reynolds_moment(1.0, x, [1.0; 3]).is_none());
            assert!(transfer(1.2, 343.0, x, [1.0; 3]).is_none());
        }
        assert!(transfer(1.2, 343.0, 127.0, [0.0; 3]).is_none());
        assert!(transfer(1.2, 343.0, 127.0, [f64::NAN; 3]).is_none());
        assert!(reynolds_moment(1.0, 1.0, [f64::INFINITY; 3]).is_none());
    }
}
