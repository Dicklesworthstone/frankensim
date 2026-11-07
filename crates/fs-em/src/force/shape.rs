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

#[cfg(test)]
mod tests {
    use super::shape_force;

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
}
