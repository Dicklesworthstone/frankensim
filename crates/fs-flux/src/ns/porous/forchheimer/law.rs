//! Isotropic quadratic pore resistance, without density or FE dependencies.
//!
//! The momentum term is beta |u-u_s| (u-u_s), beta in 1/m. Physical
//! force density multiplies this by the original fluid density. This is the
//! convention used by the material-owned directional Forchheimer calibration,
//! with no hidden factor of 1/2. An isotropic two-dimensional beta must be
//! supplied independently; a directional gas fit does not authorize it.

/// Pointwise skeleton reaction, exact velocity tangent and Newton affine load.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(super) struct Response {
    pub force: [f64; 2],
    pub tangent: [[f64; 2]; 2],
    pub load: [f64; 2],
    pub fluid_power: f64,
    pub skeleton_power: f64,
    pub dissipation: f64,
}

fn dot(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}

/// Evaluate the actual force and its C1 tangent, including the exact zero-slip
/// limit. No velocity floor, fitted permeability or clipped resistance enters.
pub(super) fn evaluate(
    beta: f64,
    velocity: [f64; 2],
    skeleton: [f64; 2],
) -> Result<Response, &'static str> {
    if !beta.is_finite() || beta < 0.0 || !velocity.into_iter().chain(skeleton).all(f64::is_finite)
    {
        return Err(
            "Forchheimer coefficient and velocities must be finite; beta must be nonnegative",
        );
    }
    if beta == 0.0 {
        return Ok(Response::default());
    }
    let relative = [velocity[0] - skeleton[0], velocity[1] - skeleton[1]];
    let speed = relative[0].hypot(relative[1]);
    if !speed.is_finite() {
        return Err("Forchheimer relative velocity is not representable");
    }
    if speed == 0.0 {
        return Ok(Response::default());
    }
    let rate = beta * speed;
    if !rate.is_finite() || rate == 0.0 {
        return Err("Forchheimer rate is not representable");
    }
    let direction = relative.map(|value| value / speed);
    let force = relative.map(|value| rate * value);
    let cross = rate * direction[0] * direction[1];
    let tangent = [
        [rate * (1.0 + direction[0] * direction[0]), cross],
        [cross, rate * (1.0 + direction[1] * direction[1])],
    ];
    // J*u - F = J*u_s + F because J*(u-u_s) = 2*F. This form
    // avoids subtracting nearly equal large terms in near co-motion.
    let load = [
        dot(tangent[0], skeleton) + force[0],
        dot(tangent[1], skeleton) + force[1],
    ];
    let response = Response {
        force,
        tangent,
        load,
        fluid_power: -dot(force, velocity),
        skeleton_power: dot(force, skeleton),
        dissipation: dot(force, relative),
    };
    if !response
        .force
        .into_iter()
        .chain(response.tangent.into_iter().flatten())
        .chain(response.load)
        .chain([
            response.fluid_power,
            response.skeleton_power,
            response.dissipation,
        ])
        .all(f64::is_finite)
        || response.dissipation <= 0.0
    {
        return Err("Forchheimer force, tangent or work is not representable");
    }
    Ok(response)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn near(actual: f64, expected: f64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance,
            "{actual:e} != {expected:e}"
        );
    }

    #[test]
    fn zero_beta_and_exact_comotion_have_no_artificial_tangent() {
        assert_eq!(
            evaluate(7.0, [0.3, -0.4], [0.3, -0.4]).unwrap(),
            Response::default()
        );
        assert_eq!(
            evaluate(0.0, [f64::MAX, 0.0], [-f64::MAX, 0.0]).unwrap(),
            Response::default()
        );
        assert_eq!(
            evaluate(0.0, [0.0; 2], [0.0; 2]).unwrap(),
            Response::default()
        );
    }

    #[test]
    fn three_four_five_force_tangent_and_signed_work() {
        let response = evaluate(2.0, [4.0, 2.0], [1.0, -2.0]).unwrap();
        assert_eq!(response.force, [30.0, 40.0]);
        near(response.tangent[0][0], 13.6, 1e-14);
        near(response.tangent[0][1], 4.8, 1e-14);
        near(response.tangent[1][0], 4.8, 1e-14);
        near(response.tangent[1][1], 16.4, 1e-14);
        near(response.load[0], 34.0, 1e-14);
        near(response.load[1], 12.0, 1e-14);
        assert_eq!(response.fluid_power, -200.0);
        assert_eq!(response.skeleton_power, -50.0);
        assert_eq!(response.dissipation, 250.0);
    }

    #[test]
    fn consistent_tangent_matches_independent_force_differences() {
        let frame = [0.6, -0.2];
        for velocity in [[0.9, 0.2], [-2.0, 0.0], [0.0, 3.0]] {
            let response = evaluate(7.3, velocity, frame).unwrap();
            assert_eq!(response.tangent[0][1], response.tangent[1][0]);
            for column in 0..2 {
                let mut plus = velocity;
                let mut minus = velocity;
                plus[column] += 1e-5;
                minus[column] -= 1e-5;
                let fp = evaluate(7.3, plus, frame).unwrap().force;
                let fm = evaluate(7.3, minus, frame).unwrap().force;
                for row in 0..2 {
                    near(
                        response.tangent[row][column],
                        (fp[row] - fm[row]) / 2e-5,
                        2e-8,
                    );
                }
            }
            for direction in [[1.0, 0.0], [0.0, 1.0], [0.4, -0.3]] {
                let product = [
                    dot(response.tangent[0], direction),
                    dot(response.tangent[1], direction),
                ];
                assert!(dot(direction, product) > 0.0);
            }
            for row in 0..2 {
                near(
                    dot(response.tangent[row], velocity) - response.load[row],
                    response.force[row],
                    1e-13,
                );
            }
        }
    }

    #[test]
    fn rotation_covariance_and_galilean_work_balance() {
        let velocity = [0.9, -0.7];
        let skeleton = [0.2, 0.4];
        let response = evaluate(3.0, velocity, skeleton).unwrap();
        let rotate = |a: [f64; 2]| [-a[1], a[0]];
        let rotated = evaluate(3.0, rotate(velocity), rotate(skeleton)).unwrap();
        assert_eq!(rotated.force, rotate(response.force));
        near(rotated.dissipation, response.dissipation, 1e-14);
        near(rotated.tangent[0][0], response.tangent[1][1], 1e-14);
        near(rotated.tangent[0][1], -response.tangent[0][1], 1e-14);
        let boost = [1.5, -0.3];
        let boosted = evaluate(
            3.0,
            [velocity[0] + boost[0], velocity[1] + boost[1]],
            [skeleton[0] + boost[0], skeleton[1] + boost[1]],
        )
        .unwrap();
        for axis in 0..2 {
            near(boosted.force[axis], response.force[axis], 1e-14);
        }
        near(boosted.dissipation, response.dissipation, 1e-14);
        near(
            boosted.fluid_power - response.fluid_power,
            -dot(response.force, boost),
            1e-14,
        );
        near(
            boosted.skeleton_power - response.skeleton_power,
            dot(response.force, boost),
            1e-14,
        );
        near(
            boosted.fluid_power + boosted.skeleton_power + boosted.dissipation,
            0.0,
            1e-14,
        );
    }

    #[test]
    fn moving_skeleton_can_power_the_fluid_without_negative_dissipation() {
        let response = evaluate(5.0, [1.0, 0.0], [2.0, 0.0]).unwrap();
        assert_eq!(response.force, [-5.0, 0.0]);
        assert_eq!(response.fluid_power, 5.0);
        assert_eq!(response.skeleton_power, -10.0);
        assert_eq!(response.dissipation, 5.0);
    }

    #[test]
    fn zero_slip_tangent_is_continuous_without_a_velocity_floor() {
        for scale in [1e-4, 1e-8, 1e-16, 1e-64] {
            let response = evaluate(2.0, [scale, 0.0], [0.0; 2]).unwrap();
            assert_eq!(response.tangent[0][0], 4.0 * scale);
            assert_eq!(response.tangent[1][1], 2.0 * scale);
            assert_eq!(response.tangent[0][1], 0.0);
            near(response.dissipation / (scale * scale * scale), 2.0, 1e-14);
        }
    }

    #[test]
    fn malformed_and_unrepresentable_laws_refuse() {
        for beta in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(evaluate(beta, [1.0, 0.0], [0.0; 2]).is_err());
        }
        assert!(evaluate(1.0, [f64::NAN, 0.0], [0.0; 2]).is_err());
        assert!(evaluate(1.0, [0.0; 2], [f64::INFINITY, 0.0]).is_err());
        assert!(evaluate(1.0, [f64::MAX, 0.0], [-f64::MAX, 0.0]).is_err());
        assert!(evaluate(f64::MAX, [2.0, 0.0], [0.0; 2]).is_err());
        assert!(evaluate(f64::MIN_POSITIVE, [1e-100, 0.0], [0.0; 2]).is_err());
        assert!(evaluate(1.0, [1e-200, 0.0], [0.0; 2]).is_err());
    }

    #[test]
    fn coefficient_and_speed_scalings_match_quadratic_force_and_cubic_loss() {
        let original = evaluate(2.0, [0.9, -0.3], [0.2, 0.1]).unwrap();
        let scaled = evaluate(6.0, [1.8, -0.6], [0.4, 0.2]).unwrap();
        for row in 0..2 {
            near(scaled.force[row], 12.0 * original.force[row], 2e-14);
            near(scaled.load[row], 12.0 * original.load[row], 2e-14);
            for column in 0..2 {
                near(
                    scaled.tangent[row][column],
                    6.0 * original.tangent[row][column],
                    2e-14,
                );
            }
        }
        near(scaled.dissipation, 24.0 * original.dissipation, 4e-14);
        near(scaled.fluid_power, 24.0 * original.fluid_power, 4e-14);
        near(scaled.skeleton_power, 24.0 * original.skeleton_power, 4e-14);
    }
}
