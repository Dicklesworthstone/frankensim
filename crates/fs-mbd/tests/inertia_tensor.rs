//! Full physical body-frame inertia, without relabelling geometry or contacts.
use fs_contact::frictional_multicontact::{FrictionalContactConfig, NormalVelocityLaw};
use fs_exec::CancelGate;
use fs_ga::Mat3;
use fs_mbd::frictional_impact::*;
use fs_mbd::{
    DynamicsError, Gravity, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState,
    UnitQuaternion, Vec3, Wrench,
};

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * a.abs().max(b.abs()).max(1.0),
        "{a:e} != {b:e}"
    );
}
fn vector(a: Vec3, b: Vec3, tolerance: f64) {
    close(a.x, b.x, tolerance);
    close(a.y, b.y, tolerance);
    close(a.z, b.z, tolerance);
}
fn basis() -> UnitQuaternion {
    UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, -0.7), 0.81).unwrap()
}
fn properties() -> MassProperties {
    MassProperties::from_principal_frame(2.0, Vec3::ZERO, Vec3::new(1.2, 1.8, 2.1), basis())
        .unwrap()
}
fn state(q: UnitQuaternion, l: Vec3) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(Vec3::new(0.3, -0.2, 0.8), q).unwrap(),
        Vec3::new(1.0, 0.2, -0.4),
        l,
    )
    .unwrap()
}
fn integrator() -> RigidBodyIntegrator {
    RigidBodyIntegrator::new(Gravity::ZERO)
}
fn world_l(s: RigidBodyState) -> Vec3 {
    s.pose()
        .orientation()
        .rotate_body_to_world(s.angular_momentum_body())
}

#[test]
fn endpoint_world_loads_close_spin_and_orbital_impulse_with_full_inertia() {
    let p = properties();
    let gravity = Gravity::new(Vec3::new(0.0, 0.0, -0.3)).unwrap();
    let integrator = RigidBodyIntegrator::new(gravity);
    let force = Vec3::new(0.7, -0.2, 0.3);
    let torque = Vec3::new(-0.1, 0.3, 0.2);
    let mut current = state(basis(), Vec3::new(0.8, -0.6, 1.1));
    let dt = 0.04;
    for _ in 0..100 {
        let wrench = |s: RigidBodyState| Wrench {
            force_world: force,
            torque_body: s.pose().orientation().rotate_world_to_body(torque),
        };
        let mut calls = 0;
        let (result, endpoint_pose) = integrator
            .step_verlet_with(current, p, wrench(current), dt, |drifted| {
                calls += 1;
                Ok::<_, DynamicsError>((wrench(drifted), drifted.pose()))
            })
            .unwrap();
        assert_eq!(calls, 1);
        assert_eq!(endpoint_pose, result.state_after.pose());
        let total_force = force.add(gravity.acceleration_world().scale(p.mass()));
        vector(
            result
                .state_after
                .linear_momentum_world()
                .sub(current.linear_momentum_world()),
            total_force.scale(dt),
            4e-14,
        );
        let midpoint_position = current
            .pose()
            .position_world()
            .scale(0.5)
            .add(result.state_after.pose().position_world().scale(0.5));
        let expected = midpoint_position.cross(total_force).add(torque).scale(dt);
        vector(
            result
                .diagnostics_after
                .angular_momentum_world
                .sub(result.diagnostics_before.angular_momentum_world),
            expected,
            2e-13,
        );
        current = result.state_after;
    }
}

#[test]
fn configuration_loaded_verlet_reverses_full_inertia_motion_without_projection() {
    let p = properties();
    let before = state(basis(), Vec3::new(0.8, -0.6, 1.1));
    // Translation spring and body-fixed dipole in a constant world field:
    // V = |x|²/2 - (R a).b. Their exact kicks depend only on configuration.
    let load = |s: RigidBodyState| Wrench {
        force_world: s.pose().position_world().scale(-1.0),
        torque_body: Vec3::new(0.4, -0.2, 0.7).cross(
            s.pose()
                .orientation()
                .rotate_world_to_body(Vec3::new(-0.3, 0.8, 0.5)),
        ),
    };
    let advance = |s| {
        integrator()
            .step_verlet_with(s, p, load(s), 0.1, |end| {
                Ok::<_, DynamicsError>((load(end), ()))
            })
            .unwrap()
            .0
            .state_after
    };
    let after = advance(before);
    let reversed = RigidBodyState::new(
        after.pose(),
        after.linear_momentum_world().scale(-1.0),
        after.angular_momentum_body().scale(-1.0),
    )
    .unwrap();
    let recovered = advance(reversed);
    vector(
        recovered.pose().position_world(),
        before.pose().position_world(),
        3e-14,
    );
    vector(
        recovered.linear_momentum_world(),
        before.linear_momentum_world().scale(-1.0),
        3e-14,
    );
    vector(
        recovered.angular_momentum_body(),
        before.angular_momentum_body().scale(-1.0),
        3e-14,
    );
    for axis in [Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0)] {
        vector(
            recovered.pose().orientation().rotate_body_to_world(axis),
            before.pose().orientation().rotate_body_to_world(axis),
            3e-14,
        );
    }
}

#[test]
fn endpoint_load_refusal_and_invalid_wrench_publish_no_candidate() {
    let before = state(basis(), Vec3::new(0.3, 0.2, -0.7));
    let saved = before;
    let failed = integrator().step_verlet_with::<(), DynamicsError>(
        before,
        properties(),
        Wrench::ZERO,
        0.1,
        |_| Err(DynamicsError::NonFinite("endpoint source")),
    );
    assert_eq!(failed, Err(DynamicsError::NonFinite("endpoint source")));
    assert!(
        integrator()
            .step_verlet_with(before, properties(), Wrench::ZERO, 0.1, |_| {
                Ok::<_, DynamicsError>((
                    Wrench {
                        force_world: Vec3::ZERO,
                        torque_body: Vec3::new(f64::NAN, 0.0, 0.0),
                    },
                    (),
                ))
            })
            .is_err()
    );
    let mut called = false;
    assert!(
        integrator()
            .step_verlet_with(before, properties(), Wrench::ZERO, 0.0, |_| {
                called = true;
                Ok::<_, DynamicsError>((Wrench::ZERO, ()))
            })
            .is_err()
    );
    assert!(!called);
    assert_eq!(before, saved);
}

#[test]
fn full_tensor_admission_retains_products_of_inertia_and_physical_triangle_law() {
    let p = properties();
    let matrix = p.inertia_tensor_body();
    assert!(matrix.m[2].abs() > 0.05 && matrix.m[5].abs() > 0.05);
    let admitted = MassProperties::from_inertia_tensor(2.0, Vec3::ZERO, matrix).unwrap();
    let w = Vec3::new(0.6, -0.4, 0.8);
    let expected = matrix.apply(fs_ga::Vec3::new(w.x, w.y, w.z));
    vector(
        admitted.angular_momentum_body_checked(w).unwrap(),
        Vec3::new(expected.x, expected.y, expected.z),
        2e-14,
    );
    vector(
        admitted
            .angular_velocity_body_checked(admitted.angular_momentum_body_checked(w).unwrap())
            .unwrap(),
        w,
        2e-14,
    );
    for (got, want) in admitted.inertia_tensor_body().m.iter().zip(matrix.m) {
        close(*got, want, 2e-14)
    }
    for invalid in [
        Mat3 {
            m: [1.0, 0.3, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
        },
        Mat3 {
            m: [1.0, 2.0, 0.0, 2.0, 1.0, 0.0, 0.0, 0.0, 1.0],
        },
        Mat3 {
            m: [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
        },
        Mat3 {
            m: [f64::NAN, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
        },
    ] {
        assert_eq!(
            MassProperties::from_inertia_tensor(2.0, Vec3::ZERO, invalid),
            Err(DynamicsError::InvalidInertiaTensor)
        );
    }
    assert_eq!(
        MassProperties::from_inertia_tensor(
            2.0,
            Vec3::ZERO,
            Mat3 {
                m: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 3.0]
            }
        ),
        Err(DynamicsError::InconsistentPrincipalInertia)
    );
    assert!(
        admitted
            .angular_momentum_body_checked(Vec3::new(f64::INFINITY, 0.0, 0.0))
            .is_err()
    );
}

#[test]
fn diagonal_admission_preserves_all_existing_step_paths_bit_for_bit() {
    let old = MassProperties::new(2.0, Vec3::ZERO, Vec3::new(1.2, 1.8, 2.1)).unwrap();
    let tensor =
        MassProperties::from_inertia_tensor(2.0, Vec3::ZERO, old.inertia_tensor_body()).unwrap();
    assert_eq!(old, tensor);
    let before = state(basis(), Vec3::new(0.3, -0.2, 0.7));
    let wrench = Wrench {
        force_world: Vec3::new(0.5, -0.2, 0.1),
        torque_body: Vec3::new(0.2, -0.1, 0.3),
    };
    assert_eq!(
        integrator().step(before, old, wrench, 0.02),
        integrator().step(before, tensor, wrench, 0.02)
    );
    assert_eq!(
        integrator().step_symplectic(before, old, wrench, 0.02),
        integrator().step_symplectic(before, tensor, wrench, 0.02)
    );
    assert_eq!(
        integrator().step_kick_drift(before, old, wrench, 0.02),
        integrator().step_kick_drift(before, tensor, wrench, 0.02)
    );
}

#[test]
fn torque_and_free_motion_are_covariant_under_a_fixed_body_coordinate_change() {
    let p = properties();
    let q = basis();
    let principal = MassProperties::new(p.mass(), Vec3::ZERO, p.principal_inertia_body()).unwrap();
    let initial = Vec3::new(0.6, -0.4, 0.8);
    for forced in [false, true] {
        for method in 0..3 {
            let mut a = state(UnitQuaternion::IDENTITY, q.rotate_body_to_world(initial));
            let mut b = state(q, initial);
            let first_world = world_l(a);
            let initial_energy = a.kinetic_energy(p).unwrap();
            let torque = if forced {
                Vec3::new(0.11, -0.07, 0.02)
            } else {
                Vec3::ZERO
            };
            for _ in 0..2000 {
                let wa = Wrench {
                    force_world: Vec3::ZERO,
                    torque_body: q.rotate_body_to_world(torque),
                };
                let wb = Wrench {
                    force_world: Vec3::ZERO,
                    torque_body: torque,
                };
                let advance = |s, properties, w| match method {
                    0 => integrator().step(s, properties, w, 0.002),
                    1 => integrator().step_symplectic(s, properties, w, 0.002),
                    _ => integrator().step_kick_drift(s, properties, w, 0.002),
                };
                a = advance(a, p, wa).unwrap().state_after;
                b = advance(b, principal, wb).unwrap().state_after;
            }
            vector(world_l(a), world_l(b), 4e-11);
            vector(
                a.pose()
                    .orientation()
                    .rotate_body_to_world(q.rotate_body_to_world(Vec3::new(0.3, 0.7, -0.1))),
                b.pose()
                    .orientation()
                    .rotate_body_to_world(Vec3::new(0.3, 0.7, -0.1)),
                4e-11,
            );
            close(
                a.kinetic_energy(p).unwrap(),
                b.kinetic_energy(principal).unwrap(),
                4e-11,
            );
            if !forced {
                vector(world_l(a), first_world, 4e-11);
                close(
                    a.kinetic_energy(p).unwrap(),
                    initial_energy,
                    if method == 0 { 4e-11 } else { 2e-7 },
                );
            }
        }
    }
}

#[test]
fn nonprincipal_symmetric_top_retains_exact_free_precession() {
    let q = basis();
    let transverse = 1.3;
    let axial = 2.2;
    let p = MassProperties::from_principal_frame(
        2.0,
        Vec3::ZERO,
        Vec3::new(transverse, transverse, axial),
        q,
    )
    .unwrap();
    let axis = q.rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0));
    let l = Vec3::new(0.8, -0.3, 1.2);
    let duration = 3.2;
    let initial = state(UnitQuaternion::IDENTITY, l);
    let result = integrator()
        .step_symplectic(initial, p, Wrench::ZERO, duration)
        .unwrap()
        .state_after;
    let space = UnitQuaternion::from_axis_angle(l, duration * l.norm_squared().sqrt() / transverse)
        .unwrap();
    let body = UnitQuaternion::from_axis_angle(
        axis,
        duration * l.dot(axis) * (1.0 / axial - 1.0 / transverse),
    )
    .unwrap();
    let test = Vec3::new(0.4, 0.8, -0.3);
    vector(
        result.pose().orientation().rotate_body_to_world(test),
        space.rotate_body_to_world(body.rotate_body_to_world(test)),
        3e-14,
    );
    vector(world_l(result), l, 2e-14);
    close(
        result.kinetic_energy(p).unwrap(),
        initial.kinetic_energy(p).unwrap(),
        2e-14,
    );
}

#[test]
fn point_impulse_uses_the_full_inverse_tensor_and_exact_midpoint_work() {
    let p = properties();
    let before = state(UnitQuaternion::IDENTITY, Vec3::new(0.4, -0.2, 0.7));
    let arm = Vec3::new(0.4, -0.2, -0.8);
    let impulse = Vec3::new(-0.3, 0.5, 0.8);
    let result = before.apply_impulse_at_body_point(p, arm, impulse).unwrap();
    vector(
        result.state_after.angular_momentum_body(),
        before.angular_momentum_body().add(arm.cross(impulse)),
        1e-14,
    );
    let expected = impulse.scale(1.0 / p.mass()).add(
        p.angular_velocity_body_checked(arm.cross(impulse))
            .unwrap()
            .cross(arm),
    );
    vector(
        result
            .point_kinematics_after
            .point_velocity_world
            .sub(result.point_kinematics_before.point_velocity_world),
        expected,
        2e-14,
    );
    close(result.work.work_energy_residual_j, 0.0, 3e-14);
    let direction = Vec3::new(0.0, 0.0, 1.0);
    let torque = arm.cross(direction);
    let effective = before
        .directional_effective_mass(p, arm, direction)
        .unwrap();
    close(
        effective.inverse_mass_kg_inverse,
        1.0 / p.mass() + torque.dot(p.angular_velocity_body_checked(torque).unwrap()),
        2e-14,
    );
}

fn controls() -> RigidImpactControls {
    RigidImpactControls {
        maximum_bodies: 4,
        maximum_contacts: 4,
        geometry_tolerance_m: 1e-10,
        velocity_tolerance_mps: 2e-8,
        linear_momentum_tolerance_ns: 1e-9,
        angular_momentum_tolerance_nms: 1e-9,
        energy_tolerance_j: 2e-8,
        contact: FrictionalContactConfig::default(),
    }
}
#[test]
fn offcentre_normal_coulomb_impacts_are_independent_of_the_geometry_frame() {
    let full = properties();
    let q = basis();
    let diagonal =
        MassProperties::new(full.mass(), Vec3::ZERO, full.principal_inertia_body()).unwrap();
    let arm = Vec3::new(0.4, -0.2, -0.8);
    for friction in [0.0, 0.35] {
        let angular = Vec3::new(0.2, -0.1, 0.3);
        let initial = RigidBodyState::new(
            Pose::new(Vec3::new(0.0, 0.0, 0.8 + 1e-12), UnitQuaternion::IDENTITY).unwrap(),
            Vec3::new(2.0, 0.8, -3.0),
            angular,
        )
        .unwrap();
        let mut a = [initial];
        let mut b = [RigidBodyState::new(
            Pose::new(initial.pose().position_world(), q).unwrap(),
            initial.linear_momentum_world(),
            q.rotate_world_to_body(angular),
        )
        .unwrap()];
        let contact = RigidFrictionalContact {
            id: 1,
            point: RigidImpactPoint {
                body: 0,
                arm_body: arm,
            },
            other: RigidImpactCounterpart::PrescribedPlane {
                point_world: Vec3::ZERO,
                velocity_world: Vec3::ZERO,
            },
            normal_world: Vec3::new(0.0, 0.0, 1.0),
            tangent_world: Vec3::new(1.0, 0.0, 0.0),
            activation_tolerance_m: 1e-10,
            friction_coefficient: friction,
            normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0.2 },
        };
        let transformed = RigidFrictionalContact {
            point: RigidImpactPoint {
                body: 0,
                arm_body: q.rotate_world_to_body(arm),
            },
            ..contact
        };
        let ra = apply_rigid_frictional_impacts(
            &mut a,
            &[full],
            &[contact],
            controls(),
            &CancelGate::new(),
        )
        .unwrap();
        let rb = apply_rigid_frictional_impacts(
            &mut b,
            &[diagonal],
            &[transformed],
            controls(),
            &CancelGate::new(),
        )
        .unwrap();
        vector(
            a[0].linear_momentum_world(),
            b[0].linear_momentum_world(),
            2e-9,
        );
        vector(world_l(a[0]), world_l(b[0]), 2e-9);
        close(ra.kinetic_after_j, rb.kinetic_after_j, 2e-9);
        close(ra.total_energy_residual_j, 0.0, 2e-9);
        if friction == 0.0 {
            let r_cross_n = arm.cross(contact.normal_world);
            let mobility = 1.0 / full.mass()
                + r_cross_n.dot(full.angular_velocity_body_checked(r_cross_n).unwrap());
            let vn = initial
                .point_kinematics(full, arm)
                .unwrap()
                .point_velocity_world
                .z;
            close(
                a[0].linear_momentum_world().z - initial.linear_momentum_world().z,
                -1.2 * vn / mobility,
                2e-9,
            );
        }
        let accepted = a;
        let gate = CancelGate::new();
        gate.request();
        assert!(
            apply_rigid_frictional_impacts(&mut a, &[full], &[contact], controls(), &gate).is_err()
        );
        assert_eq!(a, accepted);
    }
}
