//! G1/G3/G4: leading low-Re sphere lift through the original motion owner.
use super::*;
use fs_mbd::{Pose, RigidBodyState, UnitQuaternion, Vec3 as NativeVec3, Wrench};

fn state(
    body: &fs_mbd_material::principal::PrincipalBody<'_>,
    attitude: UnitQuaternion,
    velocity: NativeVec3,
    omega_body: NativeVec3,
) -> RigidBodyState {
    let inertia = body.properties().principal_inertia_body();
    RigidBodyState::new(
        Pose::new(NativeVec3::ZERO, attitude).unwrap(),
        velocity.scale(body.properties().mass()),
        NativeVec3::new(
            inertia.x * omega_body.x,
            inertia.y * omega_body.y,
            inertia.z * omega_body.z,
        ),
    )
    .unwrap()
}

fn magnitude(v: NativeVec3) -> f64 {
    v.x.hypot(v.y).hypot(v.z)
}

fn near(actual: NativeVec3, expected: NativeVec3) {
    assert!(
        magnitude(actual.sub(expected)) <= 2e-12 * magnitude(expected).max(1e-12),
        "{actual:?} != {expected:?}"
    );
}

#[test]
fn lift_uses_fluid_density_radius_cubed_and_signed_world_spin() {
    for (radius, fluid_density, solid_density) in [
        (RADIUS, 1.0, 20000.0),
        (RADIUS, 2.0, 20000.0),
        (2.0 * RADIUS, 1.0, 20000.0),
        (RADIUS, 1.0, 40000.0),
    ] {
        let source = body(solid_density, VolumeMoments::sphere(radius).unwrap());
        let body = source.principal_body(1e-12, || false).unwrap();
        let material = fluid(fluid_density, 1.0);
        let flow =
            CreepingSphere::new(&material, radius, Gravity::ZERO, NativeVec3::ZERO, limits())
                .unwrap();
        for spin in [-2.0, 0.0, 2.0] {
            let velocity = NativeVec3::new(0.1, 0.0, 0.0);
            let initial = state(
                &body,
                UnitQuaternion::IDENTITY,
                velocity,
                NativeVec3::new(0.0, 0.0, spin),
            );
            let load = flow
                .load_with_rotational_lift(&body, initial, || false)
                .unwrap();
            near(
                load.lift_force_world_n,
                NativeVec3::new(
                    0.0,
                    core::f64::consts::PI * fluid_density * radius.powi(3) * spin * 0.1,
                    0.0,
                ),
            );
            assert_eq!(load.lift_force_world_n.dot(velocity), 0.0);
            let stokes = flow.load(&body, initial, || false).unwrap();
            assert_eq!(stokes.lift_force_world_n, NativeVec3::ZERO);
            assert_eq!(load.dissipated_power_w, stokes.dissipated_power_w);
            assert_eq!(load.rotary_torque_body_nm, stokes.rotary_torque_body_nm);
        }
        let co_moving = state(
            &body,
            UnitQuaternion::IDENTITY,
            NativeVec3::ZERO,
            NativeVec3::new(0.0, 0.0, 1.0),
        );
        assert_eq!(
            flow.load_with_rotational_lift(&body, co_moving, || false)
                .unwrap()
                .lift_force_world_n,
            NativeVec3::ZERO,
        );
    }
}

#[test]
fn rotated_and_boosted_frames_preserve_lift_and_relative_power() {
    let source = sphere(20000.0);
    let body = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(1.0, 1.0);
    let relative = NativeVec3::new(0.04, -0.02, 0.01);
    let spin = NativeVec3::new(0.5, 1.0, -0.75);
    let resting =
        CreepingSphere::new(&material, RADIUS, Gravity::ZERO, NativeVec3::ZERO, limits()).unwrap();
    let reference = resting
        .load_with_rotational_lift(
            &body,
            state(&body, UnitQuaternion::IDENTITY, relative, spin),
            || false,
        )
        .unwrap();
    for attitude in [
        UnitQuaternion::IDENTITY,
        UnitQuaternion::from_axis_angle(NativeVec3::new(1.0, 2.0, -1.0), 0.7).unwrap(),
    ] {
        let ambient = NativeVec3::new(0.1, -0.2, 0.03);
        let flow =
            CreepingSphere::new(&material, RADIUS, Gravity::ZERO, ambient, limits()).unwrap();
        let slip = attitude.rotate_body_to_world(relative);
        let velocity = slip.add(ambient);
        let load = flow
            .load_with_rotational_lift(&body, state(&body, attitude, velocity, spin), || false)
            .unwrap();
        near(
            load.lift_force_world_n,
            attitude.rotate_body_to_world(reference.lift_force_world_n),
        );
        near(
            load.wrench.force_world,
            attitude.rotate_body_to_world(reference.wrench.force_world),
        );
        close(load.dissipated_power_w, reference.dissipated_power_w);
        close(
            load.background_power_w,
            load.wrench.force_world.dot(ambient),
        );
        assert!(load.lift_force_world_n.dot(slip).abs() < 1e-22);
        close(
            load.wrench.force_world.dot(velocity) + load.wrench.torque_body.dot(spin),
            load.background_power_w - load.dissipated_power_w,
        );
    }
}

#[test]
fn native_midpoint_closes_transverse_force_drive_and_reservoir_work_together() {
    let source = sphere(20000.0);
    let body = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(10.0, 1.0);
    let ambient = NativeVec3::new(0.1, -0.2, 0.03);
    let flow = CreepingSphere::new(&material, RADIUS, gravity(), ambient, limits()).unwrap();
    let attitude = UnitQuaternion::from_axis_angle(NativeVec3::new(1.0, 2.0, -1.0), 0.7).unwrap();
    let mut current = state(
        &body,
        attitude,
        ambient.add(NativeVec3::new(0.04, -0.02, 0.01)),
        NativeVec3::new(0.4, -0.8, 0.7),
    );
    let applied = Wrench {
        force_world: NativeVec3::new(0.003, -0.002, 0.001),
        torque_body: NativeVec3::new(2e-6, -3e-6, 1e-6),
    };
    let dt = 0.005;
    for _ in 0..40 {
        let result = flow
            .step_with_rotational_lift(&body, current, applied, controls(dt), || false)
            .unwrap();
        let p = body.properties();
        let after = result.dynamics.state_after;
        let mid_velocity = current
            .center_of_mass_velocity_world(p)
            .unwrap()
            .add(after.center_of_mass_velocity_world(p).unwrap())
            .scale(0.5);
        let mid_omega = p
            .angular_velocity_body_checked(current.angular_momentum_body())
            .unwrap()
            .add(
                p.angular_velocity_body_checked(after.angular_momentum_body())
                    .unwrap(),
            )
            .scale(0.5);
        let independent_lift = current
            .pose()
            .orientation()
            .rotate_body_to_world(mid_omega)
            .cross(mid_velocity.sub(ambient))
            .scale(core::f64::consts::PI * 10.0 * RADIUS.powi(3));
        near(result.load_midpoint.lift_force_world_n, independent_lift);
        let native = RigidBodyIntegrator::new(gravity())
            .step(
                current,
                p,
                Wrench {
                    force_world: result
                        .load_midpoint
                        .wrench
                        .force_world
                        .add(applied.force_world),
                    torque_body: result
                        .load_midpoint
                        .wrench
                        .torque_body
                        .add(applied.torque_body),
                },
                dt,
            )
            .unwrap();
        assert_eq!(result.dynamics, native);
        assert_eq!(
            result.load_after,
            flow.load_with_rotational_lift(&body, after, || false)
                .unwrap()
        );
        close(
            result.applied_work_j,
            dt * (applied.force_world.dot(mid_velocity) + applied.torque_body.dot(mid_omega)),
        );
        close(
            result.background_work_j,
            dt * result
                .load_midpoint
                .wrench
                .force_world
                .sub(result.load_midpoint.buoyancy_force_world_n)
                .dot(ambient),
        );
        assert!(result.viscous_heat_j > 0.0);
        assert!(result.energy_defect_j.abs() < 1e-15);
        current = after;
    }
}

#[test]
fn curved_free_trajectory_converges_to_independent_damped_spin_solution() {
    let density = 20000.0;
    let fluid_density = 50.0;
    let source = sphere(density);
    let body = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(fluid_density, 1.0);
    let flow =
        CreepingSphere::new(&material, RADIUS, Gravity::ZERO, NativeVec3::ZERO, limits()).unwrap();
    let duration = 0.4;
    let spin = 8.0;
    let speed = 0.01;
    let translation_rate = 9.0 / (2.0 * density * RADIUS.powi(2));
    let rotation_rate = 15.0 / (density * RADIUS.powi(2));
    // v_x+i*v_y = v0*exp(-lambda*t + i*theta), with the native freely
    // decaying spin omega(t)=omega0*exp(-lambda_rot*t).
    let theta = 0.75 * fluid_density / density * spin / rotation_rate
        * (-(-rotation_rate * duration).exp_m1());
    let decay = speed * (-translation_rate * duration).exp();
    let exact = NativeVec3::new(decay * theta.cos(), decay * theta.sin(), 0.0);
    let mut errors = Vec::new();
    for steps in [20, 40, 80] {
        let dt = duration / f64::from(steps);
        let mut current = state(
            &body,
            UnitQuaternion::IDENTITY,
            NativeVec3::new(speed, 0.0, 0.0),
            NativeVec3::new(0.0, 0.0, spin),
        );
        let energy = current.kinetic_energy(body.properties()).unwrap();
        let mut heat = 0.0;
        for _ in 0..steps {
            let result = flow
                .step_with_rotational_lift(&body, current, Wrench::ZERO, controls(dt), || false)
                .unwrap();
            heat += result.viscous_heat_j;
            current = result.dynamics.state_after;
        }
        let velocity = current
            .center_of_mass_velocity_world(body.properties())
            .unwrap();
        assert!(
            velocity.y > 5e-6,
            "lift must actually curve the native trajectory"
        );
        assert!((current.kinetic_energy(body.properties()).unwrap() + heat - energy).abs() < 1e-16);
        errors.push([magnitude(velocity.sub(exact)), (velocity.y - exact.y).abs()]);
    }
    for pair in errors.windows(2) {
        for component in 0..2 {
            let ratio = pair[0][component] / pair[1][component];
            assert!(
                (3.8..4.2).contains(&ratio),
                "second-order convergence: {errors:?}"
            );
        }
    }
}

#[test]
fn zero_spin_retains_stokes_motion_and_invalid_stages_publish_nothing() {
    let source = sphere(20000.0);
    let body = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(1.0, 1.0);
    let flow =
        CreepingSphere::new(&material, RADIUS, Gravity::ZERO, NativeVec3::ZERO, limits()).unwrap();
    let initial = state(
        &body,
        UnitQuaternion::IDENTITY,
        NativeVec3::new(0.01, -0.02, 0.03),
        NativeVec3::ZERO,
    );
    assert_eq!(
        flow.step(&body, initial, controls(0.005), || false)
            .unwrap(),
        flow.step_with_rotational_lift(&body, initial, Wrench::ZERO, controls(0.005), || false)
            .unwrap(),
    );
    assert_eq!(
        flow.load_with_rotational_lift(&body, initial, || true),
        Err(CreepingFlowError::Cancelled)
    );
    assert_eq!(
        flow.step_with_rotational_lift(&body, initial, Wrench::ZERO, controls(0.005), || true,),
        Err(CreepingFlowError::Cancelled)
    );
    for (velocity, omega) in [
        (NativeVec3::new(10.0, 0.0, 0.0), NativeVec3::ZERO),
        (NativeVec3::ZERO, NativeVec3::new(0.0, 0.0, 1000.0)),
    ] {
        let outside = state(&body, UnitQuaternion::IDENTITY, velocity, omega);
        assert!(matches!(
            flow.load_with_rotational_lift(&body, outside, || false),
            Err(CreepingFlowError::Reynolds { .. })
        ));
    }
    let applied = Wrench {
        force_world: NativeVec3::new(100.0, 0.0, 0.0),
        torque_body: NativeVec3::ZERO,
    };
    assert!(matches!(
        flow.step_with_rotational_lift(&body, initial, applied, controls(0.005), || false,),
        Err(CreepingFlowError::Reynolds { .. })
    ));
    assert_eq!(
        initial.linear_momentum_world(),
        NativeVec3::new(0.01, -0.02, 0.03).scale(body.properties().mass())
    );
}
