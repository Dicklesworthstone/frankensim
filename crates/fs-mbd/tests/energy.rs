//! Energy/momentum and independent Euler-top references for the rigid-body step.
use fs_mbd::{
    Gravity, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState, UnitQuaternion, Vec3,
    Wrench,
};

fn body() -> (MassProperties, RigidBodyState) {
    let properties = MassProperties::new(1.0, Vec3::ZERO, Vec3::new(1.0, 1.5, 2.0)).unwrap();
    let pose = Pose::new(
        Vec3::ZERO,
        UnitQuaternion::from_axis_angle(Vec3::new(1.0, -2.0, 3.0), 0.7).unwrap(),
    )
    .unwrap();
    let state = RigidBodyState::new(pose, Vec3::ZERO, Vec3::new(0.4, 0.7, 1.1)).unwrap();
    (properties, state)
}

fn norm(v: Vec3) -> f64 {
    v.norm_squared().sqrt()
}

#[test]
fn asymmetric_free_top_keeps_energy_and_world_momentum_over_fifty_thousand_steps() {
    let (properties, initial) = body();
    let integrator = RigidBodyIntegrator::new(Gravity::ZERO);
    let before = integrator.diagnostics(initial, properties).unwrap();
    let mut state = initial;
    let mut energy_error = 0.0_f64;
    let mut momentum_error = 0.0_f64;
    for _ in 0..50_000 {
        let step = integrator
            .step(state, properties, Wrench::ZERO, 0.05)
            .unwrap();
        state = step.state_after;
        energy_error = energy_error
            .max((step.diagnostics_after.mechanical_energy - before.mechanical_energy).abs());
        momentum_error = momentum_error.max(norm(
            step.diagnostics_after
                .angular_momentum_world
                .sub(before.angular_momentum_world),
        ));
    }
    eprintln!(
        "50,000 steps: energy drift {energy_error:e} J, world momentum drift {momentum_error:e} N m s"
    );
    assert!(
        energy_error < 2e-10 * before.mechanical_energy,
        "energy drift {energy_error:e}"
    );
    assert!(
        momentum_error < 3e-10 * norm(before.angular_momentum_world),
        "world momentum drift {momentum_error:e}"
    );
}

#[test]
fn forced_asymmetric_rotation_balances_actual_midpoint_torque_work() {
    let (properties, mut state) = body();
    let integrator = RigidBodyIntegrator::new(Gravity::new(Vec3::new(0.0, 0.0, -9.81)).unwrap());
    let initial_energy = integrator
        .diagnostics(state, properties)
        .unwrap()
        .mechanical_energy;
    let mut work = 0.0;
    for i in 0..2_000 {
        let time = f64::from(i) * 0.01;
        let wrench = Wrench {
            force_world: Vec3::new(0.03, -0.02, 9.81),
            torque_body: Vec3::new(0.07 * time.cos(), 0.05, -0.03 * time.sin()),
        };
        let step = integrator.step(state, properties, wrench, 0.01).unwrap();
        let midpoint_l = state
            .angular_momentum_body()
            .scale(0.5)
            .add(step.state_after.angular_momentum_body().scale(0.5));
        let displacement = step
            .state_after
            .pose()
            .position_world()
            .sub(state.pose().position_world());
        let supplied = wrench.force_world.dot(displacement)
            + 0.01
                * wrench
                    .torque_body
                    .dot(properties.angular_velocity_body(midpoint_l));
        let change =
            step.diagnostics_after.mechanical_energy - step.diagnostics_before.mechanical_energy;
        assert!(
            (change - supplied).abs() < 2e-13 * (1.0 + change.abs() + supplied.abs()),
            "step {i}: energy defect {:e}",
            change - supplied
        );
        work += supplied;
        state = step.state_after;
    }
    let final_energy = integrator
        .diagnostics(state, properties)
        .unwrap()
        .mechanical_energy;
    assert!((final_energy - initial_energy - work).abs() < 1e-10);
}

#[test]
fn symmetric_top_attitude_and_momentum_converge_at_second_order() {
    let properties = MassProperties::new(1.0, Vec3::ZERO, Vec3::new(2.0, 2.0, 3.5)).unwrap();
    let momentum = Vec3::new(2.2, 0.0, 2.8);
    let initial = RigidBodyState::new(Pose::identity(), Vec3::ZERO, momentum).unwrap();
    let integrator = RigidBodyIntegrator::new(Gravity::ZERO);
    let duration: f64 = 2.0;
    let phase = (3.5 - 2.0) * 0.8 / 2.0 * duration;
    let exact_l = Vec3::new(
        momentum.x * phase.cos(),
        momentum.x * phase.sin(),
        momentum.z,
    );
    // Exact symmetric-top attitude: precession about fixed world L, then the
    // body-axis correction. This checks attitude independently of the stepper.
    let exact_q = UnitQuaternion::from_axis_angle(momentum, norm(momentum) / 2.0 * duration)
        .unwrap()
        .right_exp(Vec3::new(0.0, 0.0, -phase))
        .unwrap();
    let axis = Vec3::new(1.0, 0.0, 0.0);
    let mut errors = Vec::new();
    for count in [20, 40, 80] {
        let mut state = initial;
        for _ in 0..count {
            state = integrator
                .step(state, properties, Wrench::ZERO, duration / f64::from(count))
                .unwrap()
                .state_after;
        }
        errors.push(
            norm(state.angular_momentum_body().sub(exact_l))
                + norm(
                    state
                        .pose()
                        .orientation()
                        .rotate_body_to_world(axis)
                        .sub(exact_q.rotate_body_to_world(axis)),
                ),
        );
    }
    for pair in errors.windows(2) {
        let ratio = pair[0] / pair[1];
        assert!(
            (3.5..4.5).contains(&ratio),
            "second-order ratio {ratio}, errors {errors:?}"
        );
    }
}

#[test]
fn single_axis_spin_keeps_the_exact_exponential_at_large_elapsed_time() {
    let properties = MassProperties::new(1.0, Vec3::ZERO, Vec3::new(1.0, 1.5, 2.0)).unwrap();
    let initial =
        RigidBodyState::new(Pose::identity(), Vec3::ZERO, Vec3::new(0.0, 0.0, 30.0)).unwrap();
    let step = RigidBodyIntegrator::new(Gravity::ZERO)
        .step(initial, properties, Wrench::ZERO, 25.4)
        .unwrap();
    let expected = UnitQuaternion::IDENTITY
        .right_exp(Vec3::new(0.0, 0.0, 15.0 * 25.4))
        .unwrap();
    assert_eq!(step.state_after.pose().orientation(), expected);
    assert_eq!(
        step.state_after.angular_momentum_body(),
        initial.angular_momentum_body()
    );
}

#[test]
fn asymmetric_free_rotation_is_time_reversible() {
    let (properties, initial) = body();
    let integrator = RigidBodyIntegrator::new(Gravity::ZERO);
    let forward = integrator
        .step(initial, properties, Wrench::ZERO, 0.1)
        .unwrap()
        .state_after;
    let reversed = RigidBodyState::new(
        forward.pose(),
        Vec3::ZERO,
        forward.angular_momentum_body().scale(-1.0),
    )
    .unwrap();
    let back = integrator
        .step(reversed, properties, Wrench::ZERO, 0.1)
        .unwrap()
        .state_after;
    assert!(
        norm(
            back.angular_momentum_body()
                .add(initial.angular_momentum_body())
        ) < 2e-14
    );
    let axis = Vec3::new(1.0, 0.0, 0.0);
    assert!(
        norm(
            back.pose()
                .orientation()
                .rotate_body_to_world(axis)
                .sub(initial.pose().orientation().rotate_body_to_world(axis))
        ) < 2e-14
    );
}

#[test]
fn unresolvable_rotation_refuses_without_committing_an_advance() {
    let (properties, initial) = body();
    let integrator = RigidBodyIntegrator::new(Gravity::ZERO);
    let mut state = initial;
    let result = integrator.advance(&mut state, properties, Wrench::ZERO, 1e100, 1, |_| false);
    assert_eq!(result, Err(fs_mbd::DynamicsError::RotationalSolveFailed));
    assert_eq!(state, initial);
}

#[test]
fn explicit_rotation_preserves_momentum_and_bounds_energy_in_both_lanes() {
    let (properties, initial) = body();
    let integrator = RigidBodyIntegrator::new(Gravity::ZERO);
    let reference = integrator.diagnostics(initial, properties).unwrap();
    for kick_drift in [false, true] {
        let mut state = initial;
        let mut energy_error = 0.0_f64;
        let mut momentum_error = 0.0_f64;
        for _ in 0..20_000 {
            let step = if kick_drift {
                integrator.step_kick_drift(state, properties, Wrench::ZERO, 0.05)
            } else {
                integrator.step_symplectic(state, properties, Wrench::ZERO, 0.05)
            }
            .unwrap();
            state = step.state_after;
            energy_error = energy_error.max(
                (step.diagnostics_after.mechanical_energy - reference.mechanical_energy).abs(),
            );
            momentum_error = momentum_error.max(norm(
                step.diagnostics_after
                    .angular_momentum_world
                    .sub(reference.angular_momentum_world),
            ));
        }
        eprintln!(
            "explicit kick_drift={kick_drift}: energy {energy_error:e}, momentum {momentum_error:e}"
        );
        assert!(
            momentum_error < 1e-10,
            "spatial momentum drift {momentum_error:e}"
        );
        assert!(
            energy_error < 3e-6,
            "energy left its envelope by {energy_error:e}"
        );
    }
}

#[test]
fn explicit_symmetric_top_is_exact_for_every_symmetry_axis() {
    let (_, initial) = body();
    let momentum = initial.angular_momentum_body();
    let duration = 2.3;
    let integrator = RigidBodyIntegrator::new(Gravity::ZERO);
    for axis_index in 0..3 {
        let mut inertia = [2.0; 3];
        inertia[axis_index] = 3.5;
        let properties = MassProperties::new(
            1.0,
            Vec3::ZERO,
            Vec3::new(inertia[0], inertia[1], inertia[2]),
        )
        .unwrap();
        let mut axis = [0.0; 3];
        axis[axis_index] = 1.0;
        let axis = Vec3::new(axis[0], axis[1], axis[2]);
        let axial = momentum.dot(axis);
        let phase = duration * axial * (1.0 / 3.5 - 1.0 / 2.0);
        // Independent symmetric-top formula: spatial precession then axial spin.
        let exact_q = initial
            .pose()
            .orientation()
            .right_exp(momentum.scale(duration / 2.0))
            .unwrap()
            .right_exp(axis.scale(phase))
            .unwrap();
        let parallel = axis.scale(axial);
        let exact_l = parallel
            .add(momentum.sub(parallel).scale(phase.cos()))
            .sub(axis.cross(momentum).scale(phase.sin()));
        let actual = integrator
            .step_symplectic(initial, properties, Wrench::ZERO, duration)
            .unwrap()
            .state_after;
        assert!(norm(actual.angular_momentum_body().sub(exact_l)) < 1e-12);
        for probe in [
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ] {
            assert!(
                norm(
                    actual
                        .pose()
                        .orientation()
                        .rotate_body_to_world(probe)
                        .sub(exact_q.rotate_body_to_world(probe))
                ) < 1e-12
            );
        }
    }
}

#[test]
fn explicit_body_torque_obeys_world_angular_impulse_balance() {
    let (properties, initial) = body();
    let integrator = RigidBodyIntegrator::new(Gravity::ZERO);
    let wrench = Wrench {
        force_world: Vec3::ZERO,
        torque_body: Vec3::new(0.6, -0.3, 0.8),
    };
    let h = 0.2;
    for kick_drift in [false, true] {
        let step = if kick_drift {
            integrator.step_kick_drift(initial, properties, wrench, h)
        } else {
            integrator.step_symplectic(initial, properties, wrench, h)
        }
        .unwrap();
        let torque_before = initial
            .pose()
            .orientation()
            .rotate_body_to_world(wrench.torque_body);
        let expected = if kick_drift {
            torque_before.scale(h)
        } else {
            torque_before
                .add(
                    step.state_after
                        .pose()
                        .orientation()
                        .rotate_body_to_world(wrench.torque_body),
                )
                .scale(0.5 * h)
        };
        let actual = step
            .diagnostics_after
            .angular_momentum_world
            .sub(step.diagnostics_before.angular_momentum_world);
        assert!(norm(actual.sub(expected)) < 1e-12);
    }
}

#[test]
fn explicit_asymmetric_rotation_is_reversible_and_second_order() {
    let (properties, initial) = body();
    let integrator = RigidBodyIntegrator::new(Gravity::ZERO);
    let run = |mut state, count, h| {
        for _ in 0..count {
            state = integrator
                .step_symplectic(state, properties, Wrench::ZERO, h)
                .unwrap()
                .state_after;
        }
        state
    };
    let forward = run(initial, 700, 0.03);
    let reversed = RigidBodyState::new(
        forward.pose(),
        Vec3::ZERO,
        forward.angular_momentum_body().scale(-1.0),
    )
    .unwrap();
    let back = run(reversed, 700, 0.03);
    assert!(
        norm(
            back.angular_momentum_body()
                .add(initial.angular_momentum_body())
        ) < 1e-11
    );
    let probe = Vec3::new(1.0, 0.0, 0.0);
    assert!(
        norm(
            back.pose()
                .orientation()
                .rotate_body_to_world(probe)
                .sub(initial.pose().orientation().rotate_body_to_world(probe))
        ) < 1e-11
    );

    let reference = run(initial, 4096, 2.0 / 4096.0);
    let errors: Vec<_> = [20, 40, 80]
        .into_iter()
        .map(|count| {
            let actual = run(initial, count, 2.0 / f64::from(count));
            norm(
                actual
                    .angular_momentum_body()
                    .sub(reference.angular_momentum_body()),
            ) + norm(
                actual
                    .pose()
                    .orientation()
                    .rotate_body_to_world(probe)
                    .sub(reference.pose().orientation().rotate_body_to_world(probe)),
            )
        })
        .collect();
    for pair in errors.windows(2) {
        let ratio = pair[0] / pair[1];
        assert!(
            (3.5..4.5).contains(&ratio),
            "second-order ratio {ratio}, errors {errors:?}"
        );
    }
}
