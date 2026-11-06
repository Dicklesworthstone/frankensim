//! Finite-inertia reciprocal viscous torque: recoil, heat, frames and convergence.

use fs_ga::{Mat3, So3, So3Tangent, Vec3};
use fs_time::lie::gyrostat::{HeldTorque, InertiaTensor3};
use fs_time::lie::viscous::pair::ViscousPairStepper;

fn diagonal(x: f64, y: f64, z: f64) -> Mat3 {
    Mat3 {
        m: [x, 0.0, 0.0, 0.0, y, 0.0, 0.0, 0.0, z],
    }
}

fn rotation(value: Vec3) -> So3 {
    So3::exp(So3Tangent::new(value)).unwrap()
}

fn sphere(moment: f64) -> InertiaTensor3 {
    InertiaTensor3::from_principal(Vec3::new(moment, moment, moment), So3::identity()).unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, tolerance {tolerance:e}"
    );
}

fn full_drag() -> Mat3 {
    Mat3 {
        m: [0.6, 0.1, -0.05, 0.1, 0.9, 0.2, -0.05, 0.2, 1.2],
    }
}

fn check_balances(step: &fs_time::lie::viscous::pair::ViscousPairStep, tolerance: f64) {
    assert!(step.dissipated_energy_j >= 0.0);
    close(step.work_energy_residual_j, 0.0, tolerance);
    assert!(step.angular_momentum_residual_space.norm() < tolerance);
    assert_eq!(
        step.bodies[0].coupling_angular_impulse_space
            + step.bodies[1].coupling_angular_impulse_space,
        Vec3::default(),
    );
    close(
        step.bodies[0].coupling_work_j + step.bodies[1].coupling_work_j + step.dissipated_energy_j,
        0.0,
        tolerance,
    );
    for body in step.bodies {
        close(body.work_energy_residual_j, 0.0, tolerance);
        assert!(body.angular_momentum_residual_space.norm() < tolerance);
    }
}

#[test]
fn driven_finite_pair_refines_to_analytic_recoil_phase_work_and_heat() {
    let inertias = [sphere(2.0), sphere(5.0)];
    let duration: f64 = 1.2;
    let initial = [2.1, -0.4];
    let torque = [0.7, -0.2];
    let damping: f64 = 0.8;
    let lambda = damping * (1.0 / 2.0 + 1.0 / 5.0);
    let terminal = (torque[0] / 2.0 - torque[1] / 5.0) / lambda;
    let amplitude = initial[0] - initial[1] - terminal;
    let exponential = (-lambda * duration).exp();
    let relative = terminal + amplitude * exponential;
    let relative_angle = terminal * duration + amplitude * (1.0 - exponential) / lambda;
    let total_initial = 2.0 * initial[0] + 5.0 * initial[1];
    let mean = (total_initial + (torque[0] + torque[1]) * duration) / 7.0;
    let mean_angle =
        total_initial / 7.0 * duration + 0.5 * (torque[0] + torque[1]) / 7.0 * duration * duration;
    let expected_rates = [mean + 5.0 / 7.0 * relative, mean - 2.0 / 7.0 * relative];
    let expected_angles = [
        mean_angle + 5.0 / 7.0 * relative_angle,
        mean_angle - 2.0 / 7.0 * relative_angle,
    ];
    let expected_heat = damping
        * (terminal * terminal * duration
            + 2.0 * terminal * amplitude * (1.0 - exponential) / lambda
            + amplitude * amplitude * (1.0 - exponential * exponential) / (2.0 * lambda));
    let initial_rotation = rotation(Vec3::new(0.3, -0.2, 0.4));
    let mut errors = Vec::new();
    for count in [8, 16, 32] {
        let stepper = ViscousPairStepper::new(
            inertias,
            diagonal(0.0, 0.0, damping),
            duration / count as f64,
        )
        .unwrap();
        let mut rotations = [initial_rotation; 2];
        let mut momenta = [
            Vec3::new(0.0, 0.0, 2.0 * initial[0]),
            Vec3::new(0.0, 0.0, 5.0 * initial[1]),
        ];
        let mut work = [0.0; 2];
        let mut heat = 0.0;
        let mut drift = 0.0;
        for _ in 0..count {
            let step = stepper
                .step(
                    rotations,
                    momenta,
                    torque.map(|v| HeldTorque::Body(Vec3::new(0.0, 0.0, v))),
                )
                .unwrap();
            check_balances(&step, 3e-13);
            heat += step.dissipated_energy_j;
            for i in 0..2 {
                work[i] += step.bodies[i].drive_work_j;
                drift += step.bodies[i].drift_energy_defect_j;
                rotations[i] = step.bodies[i].rotation_after;
                momenta[i] = step.bodies[i].momentum_after_body;
            }
        }
        let mut error = (heat - expected_heat).abs();
        for i in 0..2 {
            let rate = inertias[i].velocity(momenta[i]).unwrap().z;
            let expected_rotation = initial_rotation
                .compose(rotation(Vec3::new(0.0, 0.0, expected_angles[i])))
                .unwrap();
            error += (rate - expected_rates[i]).abs()
                + rotations[i]
                    .body_minus(expected_rotation)
                    .unwrap()
                    .angular
                    .norm()
                + (work[i] - torque[i] * expected_angles[i]).abs();
        }
        errors.push(error);
        close(
            momenta[0].z + momenta[1].z,
            total_initial + (torque[0] + torque[1]) * duration,
            4e-13,
        );
        let initial_energy = initial[0] * initial[0] + 2.5 * initial[1] * initial[1];
        let final_energy = 0.25 * momenta[0].z.powi(2) + 0.1 * momenta[1].z.powi(2);
        close(
            final_energy - initial_energy - work[0] - work[1] + heat - drift,
            0.0,
            8e-13,
        );
        assert!(momenta[0].z < 2.0 * initial[0] + torque[0] * duration);
        assert!(momenta[1].z > 5.0 * initial[1] + torque[1] * duration);
    }
    println!("finite-pair analytic refinement errors: {errors:?}");
    assert!(errors[1] < 0.27 * errors[0], "{errors:?}");
    assert!(errors[2] < 0.27 * errors[1], "{errors:?}");
}

type Vector = [f64; 3];
type Reference = [f64; 19]; // Two (quaternion, spin), two drive/transfer works, heat.

fn matrices() -> [Mat3; 2] {
    [
        Mat3 {
            m: [2.0, 0.2, -0.1, 0.2, 3.0, 0.15, -0.1, 0.15, 4.0],
        },
        Mat3 {
            m: [1.5, -0.18, 0.12, -0.18, 2.6, 0.2, 0.12, 0.2, 3.5],
        },
    ]
}

fn initial_rotations() -> [So3; 2] {
    [
        rotation(Vec3::new(0.2, -0.1, 0.3)),
        rotation(Vec3::new(-0.3, 0.4, -0.2)),
    ]
}

fn cross(a: Vector, b: Vector) -> Vector {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot(a: Vector, b: Vector) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn rotate(q: [f64; 4], v: Vector) -> Vector {
    let axis = [q[1], q[2], q[3]];
    let cross = cross(axis, v);
    let parallel = dot(axis, v);
    let scalar = q[0] * q[0] - dot(axis, axis);
    core::array::from_fn(|i| scalar * v[i] + 2.0 * axis[i] * parallel + 2.0 * q[0] * cross[i])
}

fn product(m: [f64; 9], v: Vector) -> Vector {
    core::array::from_fn(|i| m[3 * i] * v[0] + m[3 * i + 1] * v[1] + m[3 * i + 2] * v[2])
}

// Test-side cofactor inverse, independent of the native spectral inertia owner.
fn inverse_product(m: [f64; 9], v: Vector) -> Vector {
    let [a, b, c, _, d, e, _, _, f] = m;
    let cofactors = [
        d * f - e * e,
        c * e - b * f,
        b * e - c * d,
        c * e - b * f,
        a * f - c * c,
        b * c - a * e,
        b * e - c * d,
        b * c - a * e,
        a * d - b * b,
    ];
    let determinant = a * cofactors[0] + b * cofactors[1] + c * cofactors[2];
    product(cofactors, v).map(|x| x / determinant)
}

fn derivative(state: Reference) -> Reference {
    let q = [
        [state[0], state[1], state[2], state[3]],
        [state[7], state[8], state[9], state[10]],
    ];
    let w = [
        [state[4], state[5], state[6]],
        [state[11], state[12], state[13]],
    ];
    let inverse = q.map(|v| [v[0], -v[1], -v[2], -v[3]]);
    let second_in_first = rotate(inverse[0], rotate(q[1], w[1]));
    let relative: Vector = core::array::from_fn(|i| w[0][i] - second_in_first[i]);
    let drag = product(full_drag().m, relative);
    let coupling = [drag.map(|v| -v), rotate(inverse[1], rotate(q[0], drag))];
    let drives = [rotate(inverse[0], [0.2, -0.1, 0.3]), [-0.15, 0.22, -0.08]];
    let mut result = [0.0; 19];
    for body in 0..2 {
        let offset = 7 * body;
        let gyro = cross(w[body], product(matrices()[body].m, w[body]));
        let rhs = core::array::from_fn(|i| drives[body][i] + coupling[body][i] - gyro[i]);
        let wdot = inverse_product(matrices()[body].m, rhs);
        let imaginary = cross([q[body][1], q[body][2], q[body][3]], w[body]);
        result[offset] = -0.5 * dot([q[body][1], q[body][2], q[body][3]], w[body]);
        for i in 0..3 {
            result[offset + 1 + i] = 0.5 * (q[body][0] * w[body][i] + imaginary[i]);
            result[offset + 4 + i] = wdot[i];
        }
        result[14 + body] = dot(drives[body], w[body]);
        result[16 + body] = dot(coupling[body], w[body]);
    }
    result[18] = dot(relative, drag);
    result
}

fn reference() -> Reference {
    let mut state = [0.0; 19];
    for (i, attitude) in initial_rotations().iter().enumerate() {
        let q = attitude.as_quat();
        state[7 * i..7 * i + 4].copy_from_slice(&[q.w, q.x, q.y, q.z]);
    }
    state[4..7].copy_from_slice(&[0.6, -0.3, 0.7]);
    state[11..14].copy_from_slice(&[-0.4, 0.5, 0.2]);
    let h = 0.8 / 8192.0;
    for _ in 0..8192 {
        let k1 = derivative(state);
        let k2 = derivative(core::array::from_fn(|i| state[i] + 0.5 * h * k1[i]));
        let k3 = derivative(core::array::from_fn(|i| state[i] + 0.5 * h * k2[i]));
        let k4 = derivative(core::array::from_fn(|i| state[i] + h * k3[i]));
        state = core::array::from_fn(|i| {
            state[i] + h * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0
        });
    }
    state
}

#[test]
fn full_tensors_and_both_moving_frames_refine_against_independent_euler_rk4() {
    let inertias = matrices().map(|m| InertiaTensor3::new(m).unwrap());
    let reference = reference();
    let mut errors = Vec::new();
    for count in [8, 16, 32] {
        let stepper = ViscousPairStepper::new(inertias, full_drag(), 0.8 / count as f64).unwrap();
        let mut rotations = initial_rotations();
        let mut momenta = [
            inertias[0].momentum(Vec3::new(0.6, -0.3, 0.7)).unwrap(),
            inertias[1].momentum(Vec3::new(-0.4, 0.5, 0.2)).unwrap(),
        ];
        let mut drive_work = [0.0; 2];
        let mut coupling_work = [0.0; 2];
        let mut heat = 0.0;
        for _ in 0..count {
            let step = stepper
                .step(
                    rotations,
                    momenta,
                    [
                        HeldTorque::Space(Vec3::new(0.2, -0.1, 0.3)),
                        HeldTorque::Body(Vec3::new(-0.15, 0.22, -0.08)),
                    ],
                )
                .unwrap();
            check_balances(&step, 5e-13);
            heat += step.dissipated_energy_j;
            for i in 0..2 {
                rotations[i] = step.bodies[i].rotation_after;
                momenta[i] = step.bodies[i].momentum_after_body;
                drive_work[i] += step.bodies[i].drive_work_j;
                coupling_work[i] += step.bodies[i].coupling_work_j;
            }
        }
        let mut error = (heat - reference[18]).abs();
        for i in 0..2 {
            let offset = 7 * i;
            let q = [
                reference[offset],
                reference[offset + 1],
                reference[offset + 2],
                reference[offset + 3],
            ];
            let spin = Vec3::new(
                reference[offset + 4],
                reference[offset + 5],
                reference[offset + 6],
            );
            error += (momenta[i] - matrices()[i].apply(spin)).norm()
                + (drive_work[i] - reference[14 + i]).abs()
                + (coupling_work[i] - reference[16 + i]).abs();
            for marker in [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]] {
                let expected = rotate(q, marker);
                error += (rotations[i]
                    .rotate(Vec3::new(marker[0], marker[1], marker[2]))
                    .unwrap()
                    - Vec3::new(expected[0], expected[1], expected[2]))
                .norm();
            }
        }
        errors.push(error);
    }
    println!("finite-pair independent full-tensor RK4 refinement errors: {errors:?}");
    assert!(errors[1] < 0.27 * errors[0], "{errors:?}");
    assert!(errors[2] < 0.27 * errors[1], "{errors:?}");
}

#[test]
fn independent_material_labels_and_common_world_rotation_preserve_pair_physics() {
    let principal = [Vec3::new(2.0, 3.0, 4.0), Vec3::new(1.5, 2.5, 3.5)];
    let principal_axes = [So3::identity(), rotation(Vec3::new(0.2, 0.1, -0.3))];
    let inertias = core::array::from_fn(|i| {
        InertiaTensor3::from_principal(principal[i], principal_axes[i]).unwrap()
    });
    let basis = [
        rotation(Vec3::new(0.4, -0.3, 0.2)),
        rotation(Vec3::new(-0.2, 0.15, 0.35)),
    ];
    let changed_inertias = core::array::from_fn(|i| {
        InertiaTensor3::from_principal(principal[i], basis[i].compose(principal_axes[i]).unwrap())
            .unwrap()
    });
    let world = rotation(Vec3::new(-0.2, 0.5, 0.3));
    let changed_damping = basis[0]
        .matrix()
        .compose(full_drag())
        .compose(basis[0].matrix().transpose());
    let original = ViscousPairStepper::new(inertias, full_drag(), 0.025).unwrap();
    let changed = ViscousPairStepper::new(changed_inertias, changed_damping, 0.025).unwrap();
    let mut rotations = initial_rotations();
    let mut changed_rotations = core::array::from_fn(|i| {
        world
            .compose(rotations[i])
            .unwrap()
            .compose(basis[i].inverse())
            .unwrap()
    });
    let mut momenta = [Vec3::new(1.2, -0.9, 2.8), Vec3::new(-0.5, 1.1, 0.7)];
    let mut changed_momenta = core::array::from_fn(|i| basis[i].rotate(momenta[i]).unwrap());
    let drive = [Vec3::new(0.2, -0.1, 0.3), Vec3::new(-0.1, 0.25, 0.07)];
    for _ in 0..16 {
        let a = original
            .step(
                rotations,
                momenta,
                [HeldTorque::Space(drive[0]), HeldTorque::Body(drive[1])],
            )
            .unwrap();
        let b = changed
            .step(
                changed_rotations,
                changed_momenta,
                [
                    HeldTorque::Space(world.rotate(drive[0]).unwrap()),
                    HeldTorque::Body(basis[1].rotate(drive[1]).unwrap()),
                ],
            )
            .unwrap();
        check_balances(&a, 5e-13);
        check_balances(&b, 5e-13);
        close(a.dissipated_energy_j, b.dissipated_energy_j, 4e-13);
        for i in 0..2 {
            let aa = a.bodies[i];
            let bb = b.bodies[i];
            assert!(
                (bb.momentum_after_body - basis[i].rotate(aa.momentum_after_body).unwrap()).norm()
                    < 8e-12
            );
            let expected_rotation = world
                .compose(aa.rotation_after)
                .unwrap()
                .compose(basis[i].inverse())
                .unwrap();
            assert!(
                bb.rotation_after
                    .body_minus(expected_rotation)
                    .unwrap()
                    .angular
                    .norm()
                    < 8e-12
            );
            assert!(
                (bb.coupling_angular_impulse_space
                    - world.rotate(aa.coupling_angular_impulse_space).unwrap())
                .norm()
                    < 4e-13
            );
            close(aa.drive_work_j, bb.drive_work_j, 4e-13);
            close(aa.coupling_work_j, bb.coupling_work_j, 4e-13);
            rotations[i] = aa.rotation_after;
            changed_rotations[i] = bb.rotation_after;
            momenta[i] = aa.momentum_after_body;
            changed_momenta[i] = bb.momentum_after_body;
        }
    }
}

#[test]
fn zero_drag_preserves_both_original_held_torque_flows() {
    let inertias = matrices().map(|m| InertiaTensor3::new(m).unwrap());
    let h = 0.04;
    let pair = ViscousPairStepper::new(inertias, Mat3::zero(), h).unwrap();
    let rotations = initial_rotations();
    let momenta = [Vec3::new(1.2, -0.9, 2.8), Vec3::new(-0.5, 1.1, 0.7)];
    let drives = [
        HeldTorque::Space(Vec3::new(0.2, -0.1, 0.3)),
        HeldTorque::Body(Vec3::new(-0.1, 0.25, 0.07)),
    ];
    let actual = pair.step(rotations, momenta, drives).unwrap();
    for i in 0..2 {
        let expected = inertias[i]
            .viscous_stepper(Mat3::zero(), h)
            .unwrap()
            .step(rotations[i], momenta[i], drives[i])
            .unwrap();
        assert_eq!(actual.bodies[i].rotation_after, expected.rotation_after);
        assert_eq!(
            actual.bodies[i].momentum_after_body,
            expected.momentum_after_body
        );
        assert_eq!(
            actual.bodies[i].coupling_angular_impulse_space,
            Vec3::default()
        );
        assert_eq!(actual.bodies[i].coupling_work_j, 0.0);
        close(actual.bodies[i].drive_work_j, expected.drive_work_j, 3e-15);
        close(
            actual.bodies[i].drift_energy_defect_j,
            expected.drift_energy_defect_j,
            0.0,
        );
    }
    assert_eq!(actual.dissipated_energy_j, 0.0);
    check_balances(&actual, 5e-13);
}

#[test]
fn common_spin_has_no_fictitious_heat_and_closed_pair_recoils_without_external_work() {
    let inertias = [sphere(2.0), sphere(5.0)];
    let zero = HeldTorque::Body(Vec3::default());
    let stepper = ViscousPairStepper::new(inertias, diagonal(0.0, 0.0, 0.8), 0.05).unwrap();
    let common = stepper
        .step(
            [So3::identity(); 2],
            [Vec3::new(0.0, 0.0, 3.0), Vec3::new(0.0, 0.0, 7.5)],
            [zero; 2],
        )
        .unwrap();
    close(common.dissipated_energy_j, 0.0, 1e-28);
    assert!(common.bodies[0].coupling_angular_impulse_space.norm() < 1e-14);
    let mut rotations = [So3::identity(); 2];
    let mut momenta = [Vec3::new(0.0, 0.0, 4.0), Vec3::default()];
    let mut heat = 0.0;
    let mut drift = 0.0;
    for _ in 0..40 {
        let step = stepper.step(rotations, momenta, [zero; 2]).unwrap();
        check_balances(&step, 2e-13);
        heat += step.dissipated_energy_j;
        for i in 0..2 {
            assert_eq!(step.bodies[i].drive_work_j, 0.0);
            drift += step.bodies[i].drift_energy_defect_j;
            rotations[i] = step.bodies[i].rotation_after;
            momenta[i] = step.bodies[i].momentum_after_body;
        }
        close(momenta[0].z + momenta[1].z, 4.0, 2e-13);
        assert!(momenta[1].z > 0.0);
    }
    close(
        0.25 * momenta[0].z.powi(2) + 0.1 * momenta[1].z.powi(2) - 4.0 + heat - drift,
        0.0,
        4e-13,
    );
    assert!(heat > 0.0 && heat < 4.0);
}

#[test]
fn negative_damping_unresolved_relaxation_and_nonfinite_drives_refuse() {
    let inertias = [sphere(1.0); 2];
    for h in [0.0, -0.1, f64::INFINITY, f64::NAN] {
        assert!(ViscousPairStepper::new(inertias, full_drag(), h).is_err());
    }
    assert!(ViscousPairStepper::new(inertias, diagonal(-0.1, 0.0, 0.0), 0.1).is_err());
    assert!(
        ViscousPairStepper::new(
            inertias,
            Mat3 {
                m: [0.0, 1e-20, 0.0, 1e-20, 0.0, 0.0, 0.0, 0.0, 0.0],
            },
            0.1
        )
        .is_err()
    );
    let rotations = [So3::identity(); 2];
    let momenta = [Vec3::new(1.0, 0.0, 0.0), Vec3::default()];
    let zero = HeldTorque::Body(Vec3::default());
    let unresolved = ViscousPairStepper::new(inertias, diagonal(1.0, 0.0, 0.0), 3.0).unwrap();
    assert!(unresolved.step(rotations, momenta, [zero; 2]).is_err());
    let valid = ViscousPairStepper::new(inertias, full_drag(), 0.01).unwrap();
    assert!(
        valid
            .step(
                rotations,
                momenta,
                [zero, HeldTorque::Space(Vec3::new(f64::NAN, 0.0, 0.0))]
            )
            .is_err()
    );
    let invalid_momentum = [momenta[0], Vec3::new(0.0, f64::INFINITY, 0.0)];
    assert!(valid.step(rotations, invalid_momentum, [zero; 2]).is_err());
    check_balances(&valid.step(rotations, momenta, [zero; 2]).unwrap(), 1e-13);
}
