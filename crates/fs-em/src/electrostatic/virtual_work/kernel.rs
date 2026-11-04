//! Allocation-free cell contraction shared by native 2D and 3D field owners.
//!
//! Inputs are the original induced cell energy, physical volume, electric
//! field, displacement and one P1 shape gradient. The returned coordinate
//! force is `(U_cell I - volume E tensor D) grad(N)` [N]. Physical source and
//! isotropy admission belong to the accepted-field owners, not this kernel.

pub(crate) fn shape_force<const D: usize>(
    energy_j: f64,
    volume_m3: f64,
    electric_v_m: [f64; D],
    displacement_c_m2: [f64; D],
    gradient_per_m: [f64; D],
) -> Option<[f64; D]> {
    if !energy_j.is_finite()
        || energy_j < 0.0
        || !volume_m3.is_finite()
        || volume_m3 <= 0.0
        || electric_v_m
            .iter()
            .chain(&displacement_c_m2)
            .chain(&gradient_per_m)
            .any(|value| !value.is_finite())
    {
        return None;
    }
    let mut charge_density = 0.0;
    for axis in 0..D {
        charge_density = displacement_c_m2[axis].mul_add(gradient_per_m[axis], charge_density);
    }
    let charge = volume_m3 * charge_density;
    if !charge.is_finite() {
        return None;
    }
    let force = std::array::from_fn(|axis| {
        (-charge).mul_add(electric_v_m[axis], energy_j * gradient_per_m[axis])
    });
    force.iter().all(|value| value.is_finite()).then_some(force)
}

#[cfg(test)]
mod tests {
    use super::shape_force;

    fn near(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= 2e-8 * expected.abs().max(1e-20),
            "{actual} != {expected}"
        );
    }

    #[test]
    fn original_energy_derivative_matches_affine_stretch_and_shear_in_2d_and_3d() {
        // An affine potential on a reference simplex has this exact field.
        // Under x -> (I + t A)x it becomes (I+t A)^(-T) E; the physical
        // volume multiplies by det(I+t A). These independent closed energies
        // also include the invariant extrusion depth in the 2D case.
        let epsilon = 3e-11;
        let volume = 0.2;
        let h = 1e-5;
        for field in [[2.0_f64, -3.0, 0.0], [2.0, -3.0, 5.0]] {
            let energy = |t: f64, shear: bool| {
                if shear {
                    0.5 * epsilon
                        * volume
                        * (field[0].powi(2) + (field[1] - t * field[0]).powi(2) + field[2].powi(2))
                } else {
                    0.5 * epsilon
                        * volume
                        * (1.0 + t)
                        * ((field[0] / (1.0 + t)).powi(2) + field[1].powi(2) + field[2].powi(2))
                }
            };
            let u = energy(0.0, false);
            let d = field.map(|e| epsilon * e);
            let x_force = shape_force(u, volume, field, d, [1.0, 0.0, 0.0]).unwrap();
            let y_force = shape_force(u, volume, field, d, [0.0, 1.0, 0.0]).unwrap();
            near(
                x_force[0],
                (energy(h, false) - energy(-h, false)) / (2.0 * h),
            );
            near(y_force[0], (energy(h, true) - energy(-h, true)) / (2.0 * h));
            assert!(
                (y_force[0] - x_force[1]).abs() < 1e-25,
                "an infinitesimal rigid rotation must do zero work"
            );
            if field[2] == 0.0 {
                let planar =
                    shape_force(u, volume, [field[0], field[1]], [d[0], d[1]], [0.0, 1.0]).unwrap();
                assert_eq!(planar, [y_force[0], y_force[1]]);
            }
        }
    }

    #[test]
    fn voltage_reversal_preserves_force_and_unrepresentable_contractions_refuse() {
        let e = [2.0, -3.0, 5.0];
        let d = e.map(|value| 4e-11 * value);
        let gradient = [0.1, -0.2, 0.3];
        let a = shape_force(76e-11, 1.0, e, d, gradient).unwrap();
        let b = shape_force(76e-11, 1.0, e.map(|v| -v), d.map(|v| -v), gradient).unwrap();
        assert_eq!(a, b);
        assert!(shape_force(1.0, 1.0, [1.0; 2], [f64::MAX; 2], [2.0; 2]).is_none());
        assert!(shape_force(f64::NAN, 1.0, e, d, gradient).is_none());
    }
}
