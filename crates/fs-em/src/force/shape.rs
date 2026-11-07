//! Allocation-free contraction of an isotropic field's configurational stress.
//!
//! The original cell potential, physical volume, conjugate field vectors and
//! one P1 shape gradient give `(potential I - volume a tensor b) grad(N)` [N].
//! Electric owners use stored energy and `(E, D)`; planar magnetic A_z owners
//! use coenergy and `(H, B)`. Their distinct controls, constitutive assumptions
//! and geometry pullbacks are admitted by those owners, not by this kernel.

pub(crate) fn shape_force<const D: usize>(
    potential_j: f64,
    volume_m3: f64,
    a: [f64; D],
    b: [f64; D],
    gradient_per_m: [f64; D],
) -> Option<[f64; D]> {
    if !potential_j.is_finite()
        || potential_j < 0.0
        || !volume_m3.is_finite()
        || volume_m3 <= 0.0
        || a.iter()
            .chain(&b)
            .chain(&gradient_per_m)
            .any(|value| !value.is_finite())
    {
        return None;
    }
    let mut contraction = 0.0;
    for axis in 0..D {
        contraction = b[axis].mul_add(gradient_per_m[axis], contraction);
    }
    let integrated = volume_m3 * contraction;
    if !integrated.is_finite() {
        return None;
    }
    let force = std::array::from_fn(|axis| {
        (-integrated).mul_add(a[axis], potential_j * gradient_per_m[axis])
    });
    force.iter().all(|value| value.is_finite()).then_some(force)
}

/// A P1 motion is infinitesimally rigid iff its displacement gradient is skew.
/// Use local differences so a common translation never enters the test. The
/// allowance follows the original gradient contraction's roundoff scale;
/// there is no absolute strain deadband or user-selected physical tolerance.
pub(crate) fn rigid_planar_motion(
    displacement: [[f64; 2]; 3],
    gradients_per_m: [[f64; 2]; 3],
) -> bool {
    let mut gradient = [[0.0; 2]; 2];
    let mut scale = 0.0_f64;
    for axis in 0..2 {
        let first = displacement[1][axis] - displacement[0][axis];
        let second = displacement[2][axis] - displacement[0][axis];
        if !first.is_finite() || !second.is_finite() {
            return false;
        }
        for coordinate in 0..2 {
            let a = first * gradients_per_m[1][coordinate];
            let b = second * gradients_per_m[2][coordinate];
            let value = first.mul_add(gradients_per_m[1][coordinate], b);
            let local_scale = a.abs() + b.abs();
            if !value.is_finite() || !local_scale.is_finite() {
                return false;
            }
            gradient[axis][coordinate] = value;
            scale = scale.max(local_scale);
        }
    }
    let strain = gradient[0][0]
        .abs()
        .max(gradient[1][1].abs())
        .max((gradient[0][1] + gradient[1][0]).abs());
    strain <= (64.0 * f64::EPSILON) * scale
}

#[cfg(test)]
mod tests {
    use super::{rigid_planar_motion, shape_force};

    #[test]
    fn magnetic_pullback_differentiates_original_energy_for_linear_and_nonlinear_laws() {
        let b = [2.0_f64, -3.0];
        let volume = 0.2;
        let epsilon = 1e-5;
        for quartic in [0.0, 0.3] {
            // Independent constitutive potential w(B)=|B|^2+q|B|^4/4.
            // A_z is transported as a scalar: B -> F B / det(F), while the
            // physical volume multiplies by det(F). Mechanical work is -dW
            // at held potential coordinates, with the actual nonlinear W'.
            let density = |field: [f64; 2]| {
                let square = field[0].powi(2) + field[1].powi(2);
                square + 0.25 * quartic * square * square
            };
            let square = b[0].powi(2) + b[1].powi(2);
            let h = b.map(|value| (2.0 + quartic * square) * value);
            let coenergy = volume * (b[0] * h[0] + b[1] * h[1] - density(b));
            for shear in [false, true] {
                let energy = |t: f64| {
                    if shear {
                        volume * density([b[0] + t * b[1], b[1]])
                    } else {
                        volume * (1.0 + t) * density([b[0], b[1] / (1.0 + t)])
                    }
                };
                let gradient = if shear { [0.0, 1.0] } else { [1.0, 0.0] };
                let actual = shape_force(coenergy, volume, h, b, gradient).unwrap()[0];
                let expected = -(energy(epsilon) - energy(-epsilon)) / (2.0 * epsilon);
                assert!((actual - expected).abs() < 2e-8 * expected.abs().max(1.0));
            }
            let x = shape_force(coenergy, volume, h, b, [1.0, 0.0]).unwrap();
            let y = shape_force(coenergy, volume, h, b, [0.0, 1.0]).unwrap();
            assert!(
                (x[1] - y[0]).abs() < 2e-14,
                "rigid rotation must do no work"
            );
        }
    }

    #[test]
    fn rigid_recoil_rotation_preserves_original_storage_only_when_remanence_rotates() {
        let gradients = [[-1.0, -1.0], [1.0, 0.0], [0.0, 1.0]];
        for rate in [0.0, 1.0, -3.0, 1e-90] {
            let motion = [[0.0, 0.0], [0.0, rate], [-rate, 0.0]];
            assert!(rigid_planar_motion(motion, gradients));
        }
        assert!(rigid_planar_motion([[7.0, -9.0]; 3], gradients));
        for magnitude in [1.0, 1e-90] {
            for motion in [
                [[0.0, 0.0], [magnitude, 0.0], [0.0, magnitude]],
                [[0.0, 0.0], [0.0, 0.0], [magnitude, 0.0]],
            ] {
                assert!(!rigid_planar_motion(motion, gradients));
            }
        }
        assert!(!rigid_planar_motion([[f64::NAN, 0.0]; 3], gradients));

        let b = [0.3, -0.4];
        let br = [0.2, 0.1];
        let mu = 0.8;
        let volume = 0.2;
        let rotate = |v: [f64; 2], angle: f64| {
            let (s, c) = angle.sin_cos();
            [c * v[0] - s * v[1], s * v[0] + c * v[1]]
        };
        let energy = |angle, carried| {
            let moved_b = rotate(b, angle);
            let moved_br = if carried { rotate(br, angle) } else { br };
            0.5 * volume / mu
                * ((moved_b[0] - moved_br[0]).powi(2) + (moved_b[1] - moved_br[1]).powi(2))
        };
        let epsilon = 1e-5;
        assert!((energy(epsilon, true) - energy(-epsilon, true)).abs() < 1e-16);
        // Freezing Br in the laboratory frame introduces spurious work even
        // under a rigid rotation of this original cell and its pulled-back B.
        let frozen_work = -(energy(epsilon, false) - energy(-epsilon, false)) / (2.0 * epsilon);
        assert!((frozen_work - 0.0275).abs() < 1e-11);
    }
}
