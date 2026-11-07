//! Actual source observer consumes accepted velocities and original history.
use super::*;
use crate::tensor_material::PlaneStressPoint;
use fs_material::visco_tensor::{TensorProjectionBudget, TensorProny, TensorPronyBranch};

#[test]
fn original_material_observer_admits_complete_kelvin_history_and_actual_velocity() {
    let gate = CancelGate::new_clock_free();
    let law = TensorProny::from_tensor_shear(
        fs_material::IsotropicElastic::new(1000.0, 0.25, 0.01)
            .unwrap()
            .stiffness(),
        vec![
            TensorPronyBranch::from_gram_rows(
                core::array::from_fn(|i| core::array::from_fn(|j| if i == j { 3.0 } else { 0.0 })),
                0.01,
            )
            .unwrap(),
        ],
        0.001,
    )
    .unwrap()
    .with_kelvin_viscosity_mandel({
        let v = [0.2, -0.1, 0.3, 0.05, -0.08, 0.04];
        core::array::from_fn(|i| {
            core::array::from_fn(|j| v[i] * v[j] + if i == j { 1.0 } else { 0.0 })
        })
    })
    .unwrap();
    let image = PlaneStressPoint::new_shell(
        &law,
        &[[1.0, 0.3, 0.4, -0.2, 0.1]],
        0.03,
        TensorProjectionBudget {
            max_coordinates: 2,
            max_scalar_terms: 1_000_000,
        },
        &gate,
    )
    .unwrap()
    .viscoelastic_image(1_000_000, &gate)
    .unwrap();
    let r = image.point().memory_coordinates();
    assert_eq!(image.point().normal_strain_coordinates(), 1);
    let point = ModalPlaneStressPoint {
        modes: 1..2,
        memories: 2..2 + r,
        memory_map: None,
        strain_scale: 1.3,
        image,
    };
    let mut memory = vec![0.0; 2 + r];
    memory[..2].fill(f64::NAN); // another component's layout is not this point's history
    let q = [f64::NAN, 1e-5, f64::NAN];
    let v = [f64::NAN, 0.03, f64::NAN];
    assert!(
        point
            .image
            .normal_stress_pa(&[1.3e-5], &vec![0.0; r])
            .is_err()
    );
    point.validate(&q, &v, &memory).unwrap();
    assert!(point.validate(&q, &[0.0], &memory).is_err());
    assert!(
        point
            .validate(&q, &[0.0, f64::INFINITY, 0.0], &memory)
            .is_err()
    );
    assert!(point.validate(&[0.0, 1.0, 0.0], &v, &memory).is_err());
    memory[2] = f64::NAN;
    assert!(point.validate(&q, &v, &memory).is_err());
}
