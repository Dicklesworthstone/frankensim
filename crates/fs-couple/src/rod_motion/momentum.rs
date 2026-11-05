//! Private momentum tails across one original rod kick/drift/kick interval.
//!
//! The candidate owns the six high words. This scratch owns only their low
//! words and must be discarded before a material-loss operator changes state.
//! Rotation remains the original compensated fs-time Hamiltonian flow.

use fs_ga::{So3, Vec3};
use fs_math::dd::Dd;
use fs_time::LieStepError;

pub(super) type Vector = [Dd; 3];

#[derive(Clone, Copy, Default)]
pub(super) struct Remainder {
    pub(super) linear: [f64; 3],
    pub(super) spin: [f64; 3],
}

pub(super) fn join(high: [f64; 3], low: [f64; 3]) -> Vector {
    core::array::from_fn(|d| Dd::from_pair(high[d], low[d]))
}

pub(super) fn split(value: Vector) -> Result<([f64; 3], [f64; 3]), &'static str> {
    for component in value {
        checked(component)?;
    }
    Ok((value.map(|p| p.hi), value.map(|p| p.lo)))
}

pub(super) fn kick(momentum: Dd, force: f64, duration: f64) -> Result<Dd, &'static str> {
    checked(momentum + Dd::from_f64(duration) * Dd::from_f64(force))
}

pub(super) fn mean(before: Vector, after: Vector) -> Vector {
    core::array::from_fn(|d| before[d] * Dd::from_f64(0.5) + after[d] * Dd::from_f64(0.5))
}

// Refine the scalar quotient with its fused product residual. Generic DD
// division reconstructs denominator*quotient and may overflow at finite MAX/3.
fn divide(value: Dd, divisor: f64) -> Result<Dd, &'static str> {
    if !divisor.is_finite() || divisor <= 0.0 {
        return Err("invalid rod momentum divisor");
    }
    let quotient = value.hi / divisor;
    if !quotient.is_finite() {
        return Err("unrepresentable rod momentum quotient");
    }
    let residual = Dd::from_f64((-quotient).mul_add(divisor, value.hi)) + Dd::from_f64(value.lo);
    checked(
        Dd::from_f64(quotient)
            + Dd::from_f64(residual.hi / divisor)
            + Dd::from_f64(residual.lo / divisor),
    )
}

pub(super) fn velocity(momentum: Vector, mass: f64) -> Result<Vector, &'static str> {
    let mut result = [Dd::ZERO; 3];
    for d in 0..3 {
        result[d] = divide(momentum[d], mass)?;
    }
    Ok(result)
}

fn transform(matrix: [f64; 9], value: Vector, transpose: bool) -> Result<Vector, &'static str> {
    let mut result = [Dd::ZERO; 3];
    for i in 0..3 {
        for j in 0..3 {
            let coefficient = matrix[if transpose { 3 * j + i } else { 3 * i + j }];
            result[i] = result[i] + Dd::from_f64(coefficient) * value[j];
        }
        checked(result[i])?;
    }
    Ok(result)
}

/// Same admitted inverse inertia, applied to the retained mean spin for work.
pub(super) fn angular_velocity(
    spin: Vector,
    principal: [f64; 3],
    principal_to_body: Option<So3>,
) -> Result<Vector, &'static str> {
    let matrix = principal_to_body.map(So3::matrix).map(|m| m.m);
    let principal_spin = match matrix {
        Some(matrix) => transform(matrix, spin, true)?,
        None => spin,
    };
    let mut omega = [Dd::ZERO; 3];
    for d in 0..3 {
        omega[d] = divide(principal_spin[d], principal[d])?;
    }
    match matrix {
        Some(matrix) => transform(matrix, omega, false),
        None => Ok(omega),
    }
}

pub(super) fn work(force: [f64; 3], duration: f64, velocity: Vector) -> Result<Dd, &'static str> {
    let mut result = Dd::ZERO;
    for d in 0..3 {
        result = result + Dd::from_f64(duration) * Dd::from_f64(force[d]) * velocity[d];
    }
    checked(result)
}

/// Retain both words through cancellation against the position and COM offset.
pub(super) fn translate(
    position: f64,
    remainder: f64,
    momentum: Dd,
    duration: f64,
    mass: f64,
    centroid_offsets: Option<(f64, f64)>,
) -> Result<(f64, f64), &'static str> {
    let ratio = duration / mass;
    let by_ratio = if ratio.is_normal() {
        Some(divide(Dd::from_f64(duration), mass)? * momentum)
    } else {
        None
    };
    let drift = if let Some(drift) =
        by_ratio.filter(|value| value.hi.is_finite() && value.lo.is_finite())
    {
        drift
    } else {
        // Preserve the original velocity-first range when dt/m alone is
        // subnormal/unrepresentable or its leading product overflows before
        // a negative quotient tail can cancel (e.g. (MAX/3)*3). Refining a
        // subnormal ratio also cannot retain an underflowed quotient residual.
        divide(momentum, mass)? * Dd::from_f64(duration)
    };
    let offset =
        centroid_offsets.map_or(Dd::ZERO, |(old, new)| Dd::from_f64(old) - Dd::from_f64(new));
    // Retain the original COM-shift ordering: a large drift and opposite
    // eccentric offset may cancel before either is added to the position.
    let shift = checked(checked(drift)? + offset)?;
    let next = checked(Dd::from_pair(position, remainder) + shift)?;
    Ok((next.hi, next.lo))
}

/// Existing conservative swept-angle bound, evaluated on retained momentum.
pub(super) fn swept_angle(
    spin: Vector,
    principal: [f64; 3],
    duration: f64,
) -> Result<Dd, &'static str> {
    let scale = spin.iter().map(|p| p.hi.abs()).fold(0.0_f64, f64::max);
    if scale == 0.0 {
        return Ok(Dd::ZERO);
    }
    let mut squared = Dd::ZERO;
    for p in spin {
        let normalized = divide(p, scale)?;
        squared = squared + normalized * normalized;
    }
    let norm = checked(squared.sqrt() * Dd::from_f64(scale))?;
    checked(
        divide(norm, principal.into_iter().fold(f64::INFINITY, f64::min))? * Dd::from_f64(duration),
    )
}

pub(super) fn rotate(
    frame: So3,
    spin: Vector,
    principal: [f64; 3],
    principal_to_body: Option<So3>,
    duration: f64,
) -> Result<(So3, Vector), LieStepError> {
    fs_time::lie::compensated::rigid_body_momentum_step_in_frame(
        frame,
        spin,
        Vec3::new(principal[0], principal[1], principal[2]),
        principal_to_body.unwrap_or_else(So3::identity),
        duration,
    )
}

pub(super) fn checked(value: Dd) -> Result<Dd, &'static str> {
    if value.hi.is_finite() && value.lo.is_finite() {
        Ok(value)
    } else {
        Err("unrepresentable compensated rod momentum arithmetic")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_ga::{Quat, So3Tangent};

    const LARGE: f64 = 9_007_199_254_740_992.0;

    fn vector(x: [f64; 3]) -> Vec3 {
        Vec3::new(x[0], x[1], x[2])
    }
    fn rotation(x: [f64; 3]) -> So3 {
        So3::exp(So3Tangent::new(vector(x))).unwrap()
    }
    fn distance(a: So3, b: So3) -> f64 {
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
            .into_iter()
            .map(|p| (a.rotate(vector(p)).unwrap() - b.rotate(vector(p)).unwrap()).norm())
            .fold(0.0, f64::max)
    }

    #[test]
    fn full_tensor_work_retains_the_mean_spin_tail() {
        let before = [Dd::from_f64(LARGE), Dd::ZERO, Dd::ZERO];
        let after = [kick(before[0], 2.0, 0.5).unwrap(), Dd::ZERO, Dd::ZERO];
        // This x rotation gives the full inertia a nonzero yz entry while
        // retaining x as an exact principal axis with moment LARGE.
        let omega = angular_velocity(
            mean(before, after),
            [LARGE, 1.5 * LARGE, 2.0 * LARGE],
            Some(rotation([0.37, 0.0, 0.0])),
        )
        .unwrap();
        assert_eq!(
            work([2.0, 0.0, 0.0], 0.5, omega).unwrap(),
            Dd::from_pair(1.0, 0.5 / LARGE)
        );
    }

    fn axial_phase_oracle(axis: usize, sign: f64) -> So3 {
        // Independent angle addition: the midpoint adds one radian to LARGE.
        // This exercises the private numerical adapter; the complete rod
        // timestep separately enforces its much smaller swept-angle budget.
        let (s, c) = (0.5 * sign * LARGE).sin_cos();
        let (ds, dc) = (0.5 * sign).sin_cos();
        let mut xyz = [0.0; 3];
        xyz[axis] = s * dc + c * ds;
        So3::try_from_quat(Quat {
            w: c * dc - s * ds,
            x: xyz[0],
            y: xyz[1],
            z: xyz[2],
        })
        .unwrap()
    }

    #[test]
    fn original_rotation_adapter_uses_retained_phase_in_both_inertia_frames() {
        for (axis, principal, basis) in [
            (2, [1.0; 3], None),
            (0, [1.0, 1.5, 2.0], Some(rotation([0.37, 0.0, 0.0]))),
        ] {
            for sign in [-1.0, 1.0] {
                let mut spin = [Dd::ZERO; 3];
                spin[axis] = kick(Dd::from_f64(sign * LARGE), sign * 2.0, 0.5).unwrap();
                let (pose, retained) =
                    rotate(So3::identity(), spin, principal, basis, 1.0).unwrap();
                assert!(distance(pose, axial_phase_oracle(axis, sign)) < 5e-15);
                assert_eq!(
                    kick(retained[axis], sign * 2.0, 0.5).unwrap().to_f64(),
                    sign * (LARGE + 2.0)
                );
            }
        }
    }

    #[test]
    fn general_full_tensor_adapter_preserves_world_spin_and_reverses() {
        let basis = Some(rotation([0.3, -0.6, 0.15]));
        let initial = rotation([-0.2, 0.4, 0.1]);
        let spin = [
            Dd::from_pair(0.7, 1e-17),
            Dd::from_pair(-0.4, -1e-17),
            Dd::from_pair(0.9, 2e-17),
        ];
        let (pose, retained) = rotate(initial, spin, [2.0, 3.0, 4.0], basis, 0.3).unwrap();
        let before_world = initial.rotate(vector(spin.map(Dd::to_f64))).unwrap();
        let after_world = pose.rotate(vector(retained.map(Dd::to_f64))).unwrap();
        assert!((after_world - before_world).norm() < 2e-15);
        let (returned_pose, returned_spin) =
            rotate(pose, retained, [2.0, 3.0, 4.0], basis, -0.3).unwrap();
        assert!(distance(returned_pose, initial) < 2e-15);
        for d in 0..3 {
            assert!((returned_spin[d] - spin[d]).to_f64().abs() < 5e-29);
        }
    }

    #[test]
    fn rotation_budget_includes_a_midpoint_tail_hidden_from_the_high_word() {
        let spin = [Dd::ZERO, Dd::ZERO, Dd::from_pair(LARGE, 1.0)];
        let swept = swept_angle(spin, [1.0; 3], 1.0 / LARGE).unwrap();
        assert_eq!(swept, Dd::from_pair(1.0, 1.0 / LARGE));
        assert!(Dd::ONE.lt(swept));
    }

    #[test]
    fn original_finite_translation_range_survives_quotient_refinement() {
        let maximum = translate(0.0, 0.0, Dd::from_f64(3.0), f64::MAX, 3.0, None).unwrap();
        assert_eq!(maximum, (f64::MAX, 0.0));
        let tiny = f64::from_bits(2);
        assert_eq!(
            translate(0.0, 0.0, Dd::from_f64(3.0), tiny, 3.0, None).unwrap(),
            (tiny, 0.0)
        );
        assert_eq!(
            translate(
                f64::MAX,
                0.0,
                Dd::from_f64(f64::MAX),
                1.0,
                1.0,
                Some((0.0, f64::MAX))
            )
            .unwrap(),
            (f64::MAX, 0.0)
        );
    }

    #[test]
    fn quotient_tail_and_eccentric_offset_survive_position_cancellation() {
        let q = f64::MAX / 3.0;
        let expected = (-q).mul_add(3.0, f64::MAX) / 3.0;
        let (actual, _) = translate(-q, 0.0, Dd::ONE, f64::MAX, 3.0, None).unwrap();
        assert!((actual - expected).abs() <= 8.0 * f64::EPSILON * expected.abs());
        let large = 9_007_199_254_740_992.0;
        assert_eq!(
            translate(
                -large,
                0.125,
                Dd::from_pair(large, 1.0),
                1.0,
                1.0,
                Some((1.0, -1.0))
            )
            .unwrap(),
            (3.125, 0.0)
        );
    }

    #[test]
    fn retained_kicks_and_mean_work_share_the_same_momentum() {
        let large = 9_007_199_254_740_992.0;
        for sign in [-1.0, 1.0] {
            let before = [Dd::from_f64(sign * large), Dd::ZERO, Dd::ZERO];
            let middle = [
                kick(before[0], 2.0 * sign, 0.5).unwrap(),
                Dd::ZERO,
                Dd::ZERO,
            ];
            let after = [
                kick(middle[0], 2.0 * sign, 0.5).unwrap(),
                Dd::ZERO,
                Dd::ZERO,
            ];
            assert_eq!(after[0].to_f64(), sign * (large + 2.0));
            let first = work(
                [2.0 * sign, 0.0, 0.0],
                0.5,
                velocity(mean(before, middle), 1.0).unwrap(),
            )
            .unwrap();
            let second = work(
                [2.0 * sign, 0.0, 0.0],
                0.5,
                velocity(mean(middle, after), 1.0).unwrap(),
            )
            .unwrap();
            assert_eq!(first + second, Dd::from_pair(2.0 * large, 2.0));
        }
    }
}
