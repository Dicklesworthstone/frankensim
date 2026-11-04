//! G1/G3: rigid ambient rotation, local COM transport, work and frame covariance.
use super::*;
use fs_mbd::{MassProperties, RigidBodyState};
use fs_mbd_material::ellipsoid_flow::RigidAmbientFlow;

fn set_velocities(
    state: RigidBodyState,
    properties: MassProperties,
    velocity_world: Vec3,
    omega_world: Vec3,
) -> RigidBodyState {
    let omega = state.pose().orientation().rotate_world_to_body(omega_world);
    let inertia = properties.principal_inertia_body();
    RigidBodyState::new(
        state.pose(),
        velocity_world.scale(properties.mass()),
        Vec3::new(
            inertia.x * omega.x,
            inertia.y * omega.y,
            inertia.z * omega.z,
        ),
    )
    .unwrap()
}

#[test]
fn local_translation_and_relative_spin_use_actual_com_and_authored_shape_axes() {
    let boundary = EllipsoidBoundary {
        semiaxes_m: GaVec3::new(0.01, 0.025, 0.016),
        body_from_ellipsoid: Se3::from_parts(
            rotation(0.2, -0.3, 0.1),
            GaVec3::new(0.04, -0.02, 0.03),
        )
        .unwrap(),
    };
    let source = body(boundary, 20000.0);
    let native = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(1.0, 1.0);
    let reference = Vec3::new(0.2, -0.1, 0.3);
    let reference_velocity = Vec3::new(0.01, -0.02, 0.005);
    let omega = Vec3::new(0.12, -0.08, 0.05);
    let ambient = RigidAmbientFlow::new(reference, reference_velocity, omega).unwrap();
    let medium = CreepingEllipsoid::new(
        &material,
        boundary,
        Gravity::ZERO,
        Vec3::ZERO,
        limits(),
        quadrature(),
        || false,
    )
    .unwrap()
    .with_rigid_ambient(ambient)
    .unwrap();
    assert_eq!(medium.ambient_flow(), ambient);
    let pose = Se3::from_parts(rotation(-0.4, 0.25, 0.15), GaVec3::new(0.3, 0.2, -0.1)).unwrap();
    let world_from_shape = pose.compose(boundary.body_from_ellipsoid).unwrap();
    let center = mv(world_from_shape.translation());
    let seed = native.state_from_kinematics(pose, Twist::zero()).unwrap();
    vector_close(seed.pose().position_world(), center, 2e-16);
    // Explicit source field at the actual offset COM, not the authored origin.
    let flow_velocity = reference_velocity.add(omega.cross(center.sub(reference)));
    let slip = Vec3::new(0.025, -0.013, 0.018);
    let spin_slip = Vec3::new(-0.04, 0.07, -0.03);
    let state = set_velocities(
        seed,
        native.properties(),
        flow_velocity.add(slip),
        omega.add(spin_slip),
    );
    let load = medium.load(&native, state, || false).unwrap();
    let rotation = world_from_shape.rotation();
    let slip_axes = rotation.inverse().rotate(gv(slip)).unwrap();
    let spin_axes = rotation.inverse().rotate(gv(spin_slip)).unwrap();
    let c = medium.resistance().translation_kg_s;
    let d = medium.resistance().rotation_kg_m2_s;
    let force = mv(rotation
        .rotate(GaVec3::new(
            -c[0] * slip_axes.x,
            -c[1] * slip_axes.y,
            -c[2] * slip_axes.z,
        ))
        .unwrap());
    let torque_world = mv(rotation
        .rotate(GaVec3::new(
            -d[0] * spin_axes.x,
            -d[1] * spin_axes.y,
            -d[2] * spin_axes.z,
        ))
        .unwrap());
    let torque_body = state
        .pose()
        .orientation()
        .rotate_world_to_body(torque_world);
    vector_close(load.drag_force_world_n, force, 3e-16);
    vector_close(load.rotary_torque_body_nm, torque_body, 3e-18);
    relative(
        load.dissipated_power_w,
        -force.dot(slip) - torque_world.dot(spin_slip),
        2e-13,
    );
    relative(
        load.background_power_w,
        force.dot(flow_velocity) + torque_world.dot(omega),
        2e-13,
    );
    let actual_omega = native
        .properties()
        .angular_velocity_body_checked(state.angular_momentum_body())
        .unwrap();
    let body_power = load.drag_force_world_n.dot(
        state
            .center_of_mass_velocity_world(native.properties())
            .unwrap(),
    ) + load.rotary_torque_body_nm.dot(actual_omega);
    assert!((body_power + load.dissipated_power_w - load.background_power_w).abs() < 1e-17);

    let co_moving = set_velocities(seed, native.properties(), flow_velocity, omega);
    let zero = medium.load(&native, co_moving, || false).unwrap();
    assert!(magnitude(zero.drag_force_world_n) < 2e-17);
    assert!(magnitude(zero.rotary_torque_body_nm) < 2e-20);
    assert!(zero.dissipated_power_w < 1e-30);
    assert!(zero.rotation_reynolds > 0.0);

    // A different reference point describes the same field if its velocity is moved too.
    let shift = Vec3::new(-0.17, 0.08, 0.09);
    let equivalent = RigidAmbientFlow::new(
        reference.add(shift),
        reference_velocity.add(omega.cross(shift)),
        omega,
    )
    .unwrap();
    let second = medium
        .with_rigid_ambient(equivalent)
        .unwrap()
        .load(&native, state, || false)
        .unwrap();
    vector_close(load.drag_force_world_n, second.drag_force_world_n, 3e-17);
    vector_close(
        load.rotary_torque_body_nm,
        second.rotary_torque_body_nm,
        3e-20,
    );
    relative(load.background_power_w, second.background_power_w, 2e-13);
}

// Independent converged exponential-series evaluation of the constant 4x4
// sphere equation x'=v, v'=-lambda*(v-Omega cross x). There is no midpoint solve.
fn exact_planar_sphere(initial: [f64; 4], lambda: f64, omega: f64, time: f64) -> [f64; 4] {
    let mut sum = initial;
    let mut term = initial;
    for n in 1..=48 {
        let product = [
            term[2],
            term[3],
            -lambda * (term[2] + omega * term[1]),
            -lambda * (term[3] - omega * term[0]),
        ];
        for k in 0..4 {
            term[k] = product[k] * (time / n as f64);
            sum[k] += term[k];
        }
    }
    sum
}

#[test]
fn native_rotating_flow_trajectory_converges_to_independent_sphere_translation_and_spin() {
    let radius = 0.01;
    let boundary = shape(radius, radius, radius);
    let source = body(boundary, 20000.0);
    let native = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(1.0, 1.0);
    let omega = 0.2;
    let medium = CreepingEllipsoid::new(
        &material,
        boundary,
        Gravity::ZERO,
        Vec3::ZERO,
        limits(),
        quadrature(),
        || false,
    )
    .unwrap()
    .with_rigid_ambient(
        RigidAmbientFlow::new(Vec3::ZERO, Vec3::ZERO, Vec3::new(0.0, 0.0, omega)).unwrap(),
    )
    .unwrap();
    let initial = [0.04, -0.03, 0.01, 0.02];
    let seed = native
        .state_from_kinematics(
            Se3::from_parts(So3::identity(), GaVec3::new(initial[0], initial[1], 0.0)).unwrap(),
            Twist::zero(),
        )
        .unwrap();
    let seed = set_velocities(
        seed,
        native.properties(),
        Vec3::new(initial[2], initial[3], 0.0),
        Vec3::ZERO,
    );
    let time = 0.4;
    // Independent spherical Stokes values, not the prepared coefficient report.
    let lambda = 6.0 * core::f64::consts::PI * radius / native.properties().mass();
    let rotary_rate = 8.0 * core::f64::consts::PI * radius.powi(3)
        / native.properties().principal_inertia_body().z;
    let exact = exact_planar_sphere(initial, lambda, omega, time);
    let exact_spin = omega * (1.0 - (-rotary_rate * time).exp());
    let mut errors = Vec::new();
    for steps in [10, 20, 40] {
        let dt = time / steps as f64;
        let mut state = seed;
        let mut balance = 0.0;
        let mut background = 0.0;
        let mut heat = 0.0;
        for _ in 0..steps {
            let step = medium.step(&native, state, controls(dt), || false).unwrap();
            let flow = step.flow;
            let midpoint = state.pose().position_world().add(
                flow.dynamics
                    .state_after
                    .pose()
                    .position_world()
                    .sub(state.pose().position_world())
                    .scale(0.5),
            );
            let velocity_mid = state
                .center_of_mass_velocity_world(native.properties())
                .unwrap()
                .scale(0.5)
                .add(
                    flow.dynamics
                        .state_after
                        .center_of_mass_velocity_world(native.properties())
                        .unwrap()
                        .scale(0.5),
                );
            let expected_force = velocity_mid
                .sub(Vec3::new(-omega * midpoint.y, omega * midpoint.x, 0.0))
                .scale(-6.0 * core::f64::consts::PI * radius);
            vector_close(flow.load_midpoint.drag_force_world_n, expected_force, 1e-16);
            balance += flow.stored_energy_change_j + flow.viscous_heat_j - flow.background_work_j;
            background += flow.background_work_j;
            heat += flow.viscous_heat_j;
            state = flow.dynamics.state_after;
            assert_eq!(step.ambient_rotation_increment_rad, dt * omega);
        }
        let velocity = state
            .center_of_mass_velocity_world(native.properties())
            .unwrap();
        let spin = native
            .properties()
            .angular_velocity_body_checked(state.angular_momentum_body())
            .unwrap();
        let position = state.pose().position_world();
        errors.push([
            (position.x - exact[0]).hypot(position.y - exact[1]),
            (velocity.x - exact[2]).hypot(velocity.y - exact[3]),
            (spin.z - exact_spin).abs(),
        ]);
        assert!(balance.abs() < 1e-15);
        assert!(heat > 0.0 && background != 0.0);
    }
    for pair in errors.windows(2) {
        for axis in 0..3 {
            assert!(
                pair[0][axis] > 3.7 * pair[1][axis],
                "refinement: {errors:?}"
            );
        }
    }
}

#[test]
fn triaxial_native_flow_is_covariant_and_closes_rotating_reservoir_work() {
    let boundary = EllipsoidBoundary {
        semiaxes_m: GaVec3::new(0.01, 0.025, 0.016),
        body_from_ellipsoid: Se3::from_parts(
            rotation(0.2, -0.3, 0.1),
            GaVec3::new(0.04, -0.02, 0.03),
        )
        .unwrap(),
    };
    let source = body(boundary, 20000.0);
    let native = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(1.0, 1.0);
    let rotate = rotation(-0.5, 0.4, 0.3);
    let shift = Vec3::new(0.7, -0.4, 0.6);
    let map = |v| mv(rotate.rotate(gv(v)).unwrap());
    let change = Se3::from_parts(rotate, gv(shift)).unwrap();
    let origin = Vec3::new(0.2, -0.1, 0.3);
    let velocity = Vec3::new(0.01, -0.02, 0.005);
    let omega = Vec3::new(0.12, -0.08, 0.05);
    let gravity = Gravity::new(Vec3::new(0.001, -0.002, 0.0)).unwrap();
    let first = CreepingEllipsoid::new(
        &material,
        boundary,
        gravity,
        Vec3::ZERO,
        limits(),
        quadrature(),
        || false,
    )
    .unwrap()
    .with_rigid_ambient(RigidAmbientFlow::new(origin, velocity, omega).unwrap())
    .unwrap();
    let second = CreepingEllipsoid::new(
        &material,
        boundary,
        Gravity::new(map(gravity.acceleration_world())).unwrap(),
        Vec3::ZERO,
        limits(),
        quadrature(),
        || false,
    )
    .unwrap()
    .with_rigid_ambient(
        RigidAmbientFlow::new(map(origin).add(shift), map(velocity), map(omega)).unwrap(),
    )
    .unwrap();
    let pose = Se3::from_parts(rotation(0.1, 0.2, -0.3), GaVec3::new(0.3, 0.2, 0.1)).unwrap();
    let twist = Twist::new(
        GaVec3::new(0.12, -0.07, 0.09),
        GaVec3::new(0.03, 0.02, -0.01),
    );
    let mut state = native.state_from_kinematics(pose, twist).unwrap();
    let mut rotated = native
        .state_from_kinematics(change.compose(pose).unwrap(), twist)
        .unwrap();
    let applied = Wrench {
        force_world: Vec3::new(0.004, -0.003, 0.002),
        torque_body: Vec3::new(2e-5, -3e-5, 1e-5),
    };
    let mut totals = [0.0; 4];
    for _ in 0..32 {
        let step = first
            .step_with_wrench(&native, state, applied, controls(0.005), || false)
            .unwrap();
        let rotated_step = second
            .step_with_wrench(
                &native,
                rotated,
                Wrench {
                    force_world: map(applied.force_world),
                    ..applied
                },
                controls(0.005),
                || false,
            )
            .unwrap();
        let a = step.flow;
        let b = rotated_step.flow;
        state = a.dynamics.state_after;
        rotated = b.dynamics.state_after;
        vector_close(
            rotated.linear_momentum_world(),
            map(state.linear_momentum_world()),
            2e-13,
        );
        vector_close(
            rotated.pose().position_world(),
            map(state.pose().position_world()).add(shift),
            2e-13,
        );
        vector_close(
            rotated.angular_momentum_body(),
            state.angular_momentum_body(),
            2e-15,
        );
        relative(a.viscous_heat_j, b.viscous_heat_j, 2e-11);
        assert!((a.background_work_j - b.background_work_j).abs() < 2e-16);
        totals[0] += a.stored_energy_change_j;
        totals[1] += a.viscous_heat_j;
        totals[2] += a.background_work_j;
        totals[3] += a.applied_work_j;
    }
    assert!(totals[1] > 0.0);
    assert!((totals[0] + totals[1] - totals[2] - totals[3]).abs() < 2e-12);
}

#[test]
fn uniform_flow_keeps_its_native_path_and_rotating_flow_refuses_atomically() {
    let boundary = shape(0.01, 0.02, 0.015);
    let source = body(boundary, 20000.0);
    let native = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(1.0, 1.0);
    let velocity = Vec3::new(0.01, -0.02, 0.015);
    let prepare = || {
        CreepingEllipsoid::new(
            &material,
            boundary,
            Gravity::ZERO,
            velocity,
            limits(),
            quadrature(),
            || false,
        )
        .unwrap()
    };
    let first = prepare();
    let uniform = prepare()
        .with_rigid_ambient(
            RigidAmbientFlow::new(Vec3::new(1e100, -1e100, 1e100), velocity, Vec3::ZERO).unwrap(),
        )
        .unwrap();
    let state = native
        .state_from_kinematics(
            Se3::identity(),
            Twist::new(
                GaVec3::new(0.1, -0.07, 0.05),
                GaVec3::new(0.03, 0.02, -0.01),
            ),
        )
        .unwrap();
    let applied = Wrench {
        force_world: Vec3::new(0.002, -0.001, 0.003),
        torque_body: Vec3::new(1e-5, -2e-5, 3e-5),
    };
    assert_eq!(
        first
            .step_with_wrench(&native, state, applied, controls(0.005), || false)
            .unwrap(),
        uniform
            .step_with_wrench(&native, state, applied, controls(0.005), || false)
            .unwrap(),
    );
    let flowing = prepare()
        .with_rigid_ambient(
            RigidAmbientFlow::new(Vec3::ZERO, velocity, Vec3::new(0.3, -0.2, 0.1)).unwrap(),
        )
        .unwrap();
    let mut polls = 0;
    let accepted = flowing
        .step(&native, state, controls(0.005), || {
            polls += 1;
            false
        })
        .unwrap();
    let mut at = 0;
    assert_eq!(
        flowing.step(&native, state, controls(0.005), || {
            at += 1;
            at == polls
        }),
        Err(EllipsoidFlowError::Cancelled)
    );
    assert_eq!(
        flowing
            .step(&native, state, controls(0.005), || false)
            .unwrap(),
        accepted
    );
    let mut angle = controls(0.005);
    angle.maximum_rotation_increment_rad = 1e-4;
    assert!(flowing.step(&native, state, angle, || false).is_err());
    let mut budget = controls(0.005);
    budget.maximum_coupling_iterations = 1;
    assert!(flowing.step(&native, state, budget, || false).is_err());
    assert!(
        prepare()
            .with_rigid_ambient(
                RigidAmbientFlow::new(Vec3::ZERO, Vec3::ZERO, Vec3::new(1000.0, 0.0, 0.0)).unwrap(),
            )
            .is_err()
    );
    assert!(RigidAmbientFlow::new(Vec3::ZERO, Vec3::ZERO, Vec3::new(f64::NAN, 0.0, 0.0)).is_err());
    assert!(
        RigidAmbientFlow::new(Vec3::new(f64::INFINITY, 0.0, 0.0), Vec3::ZERO, Vec3::ZERO).is_err()
    );
    assert!(
        flowing
            .ambient_flow()
            .velocity_at(Vec3::new(f64::NAN, 0.0, 0.0))
            .is_err()
    );
}
