//! Original material and interface data drive exact native velocity relaxation.
use super::*;
use fs_mbd::{Pose, RigidBodyState, Vec3 as NativeVec3, Wrench};
use fs_mbd_material::creeping_flow::exact::ExactCreepingStepControls;
use fs_mbd_material::creeping_flow::navier_slip::NavierSlipBoundary;

fn exact_controls(duration_s: f64) -> ExactCreepingStepControls {
    ExactCreepingStepControls {
        duration_s,
        maximum_energy_defect_j: 1e-11,
    }
}
fn initial(
    body: &fs_mbd_material::principal::PrincipalBody<'_>,
    velocity: NativeVec3,
    omega: NativeVec3,
) -> RigidBodyState {
    RigidBodyState::new(
        Pose::identity(),
        velocity.scale(body.properties().mass()),
        body.properties()
            .angular_momentum_body_checked(omega)
            .unwrap(),
    )
    .unwrap()
}
fn close_vec(actual: NativeVec3, expected: NativeVec3) {
    close(actual.x, expected.x);
    close(actual.y, expected.y);
    close(actual.z, expected.z);
}

#[test]
fn material_swaps_reach_exact_settling_and_spin_beyond_midpoint_step_limits() {
    let ambient = NativeVec3::new(0.03, -0.02, 0.01);
    let velocity = NativeVec3::new(0.05, 0.01, 0.04);
    let omega = NativeVec3::new(0.3, -0.4, 0.5);
    for (density, fluid_density, viscosity) in [
        (10000.0, 1.0, 1.0),
        (20000.0, 1.0, 1.0),
        (20000.0, 2.0, 2.0),
    ] {
        let source = sphere(density);
        let body = source.principal_body(1e-12, || false).unwrap();
        let material = fluid(fluid_density, viscosity);
        let reservoir =
            CreepingSphere::new(&material, RADIUS, gravity(), ambient, limits()).unwrap();
        let state = initial(&body, velocity, omega);
        let duration = 5.0;
        assert!(
            reservoir
                .step(&body, state, controls(duration), || false)
                .is_err()
        );
        let step = reservoir
            .step_exact(&body, state, exact_controls(duration), || false)
            .unwrap();
        let mass = body.properties().mass();
        let volume = 4.0 * core::f64::consts::PI * RADIUS.powi(3) / 3.0;
        let c = 6.0 * core::f64::consts::PI * viscosity * RADIUS;
        let c_rot = 8.0 * core::f64::consts::PI * viscosity * RADIUS.powi(3);
        let tau = mass / c;
        let tau_rot = body.properties().principal_inertia_body().x / c_rot;
        let terminal = ambient.add(
            gravity()
                .acceleration_world()
                .scale((mass - fluid_density * volume) / c),
        );
        let e = (-duration / tau).exp();
        let expected_velocity = terminal.add(velocity.sub(terminal).scale(e));
        let expected_position = terminal.scale(duration).add(
            velocity
                .sub(terminal)
                .scale(tau * -(-duration / tau).exp_m1()),
        );
        close_vec(
            step.dynamics
                .state_after
                .center_of_mass_velocity_world(body.properties())
                .unwrap(),
            expected_velocity,
        );
        close_vec(
            step.dynamics.state_after.pose().position_world(),
            expected_position,
        );
        close_vec(
            body.properties()
                .angular_velocity_body_checked(step.dynamics.state_after.angular_momentum_body())
                .unwrap(),
            omega.scale((-duration / tau_rot).exp()),
        );
        assert!(step.viscous_heat_j > 0.0);
        assert!(step.energy_defect_j.abs() < 1e-11);
        assert!(step.diffusion_ratio <= limits().maximum_diffusion_ratio);
        assert!(core::ptr::eq(reservoir.source(), &material));
        close(
            step.stored_energy_change_j + step.viscous_heat_j,
            step.background_work_j,
        );
    }
}

#[test]
fn navier_source_changes_the_exact_translation_and_rotation_resistances() {
    let source = sphere(20000.0);
    let body = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(1.0, 1.0);
    let ambient = NativeVec3::new(0.03, -0.01, 0.02);
    let velocity = NativeVec3::new(0.1, -0.2, 0.05);
    let omega = NativeVec3::new(0.0, 0.0, 1.0);
    let state = initial(&body, velocity, omega);
    for (slip, translation_factor, rotation_factor) in [(0.0, 1.0, 1.0), (RADIUS, 0.75, 0.25)] {
        let interface = super::navier_slip_source::constant(
            source.regions()[0].material().material(),
            slip,
            400.0,
        );
        let reservoir = CreepingSphere::new(&material, RADIUS, Gravity::ZERO, ambient, limits())
            .unwrap()
            .with_navier_slip(NavierSlipBoundary::new(&interface).unwrap())
            .unwrap();
        let duration = 2.0;
        let step = reservoir
            .step_exact(&body, state, exact_controls(duration), || false)
            .unwrap();
        let rate =
            translation_factor * 6.0 * core::f64::consts::PI * RADIUS / body.properties().mass();
        let spin_rate = rotation_factor * 8.0 * core::f64::consts::PI * RADIUS.powi(3)
            / body.properties().principal_inertia_body().x;
        close_vec(
            step.dynamics
                .state_after
                .center_of_mass_velocity_world(body.properties())
                .unwrap(),
            ambient.add(velocity.sub(ambient).scale((-rate * duration).exp())),
        );
        close_vec(
            body.properties()
                .angular_velocity_body_checked(step.dynamics.state_after.angular_momentum_body())
                .unwrap(),
            omega.scale((-spin_rate * duration).exp()),
        );
        let expected_heat = 0.5
            * body.properties().mass()
            * velocity.sub(ambient).norm_squared()
            * -(-2.0 * rate * duration).exp_m1()
            + 0.5
                * body.properties().principal_inertia_body().x
                * -(-2.0 * spin_rate * duration).exp_m1();
        close(step.viscous_heat_j, expected_heat);
        assert!(core::ptr::eq(
            reservoir.navier_slip().unwrap().source(),
            &interface
        ));
    }
}

#[test]
fn held_drive_reaches_its_terminal_slip_and_separates_source_work_from_heat() {
    let source = sphere(20000.0);
    let body = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(1.0, 1.0);
    let ambient = NativeVec3::new(0.03, -0.02, 0.01);
    let terminal_slip = NativeVec3::new(0.1, -0.2, 0.03);
    let c = 6.0 * core::f64::consts::PI * RADIUS;
    let applied = Wrench {
        force_world: terminal_slip.scale(c),
        torque_body: NativeVec3::ZERO,
    };
    let reservoir =
        CreepingSphere::new(&material, RADIUS, Gravity::ZERO, ambient, limits()).unwrap();
    let state = initial(&body, ambient, NativeVec3::ZERO);
    let duration = 8.0;
    let step = reservoir
        .step_exact_with_wrench(&body, state, applied, exact_controls(duration), || false)
        .unwrap();
    let rate = c / body.properties().mass();
    close_vec(
        step.dynamics
            .state_after
            .center_of_mass_velocity_world(body.properties())
            .unwrap(),
        ambient.add(terminal_slip.scale(-(-rate * duration).exp_m1())),
    );
    close(
        step.applied_work_j,
        applied.force_world.dot(step.dynamics.displacement_world),
    );
    close(
        step.background_work_j,
        step.dynamics.drag_impulse_world.dot(ambient),
    );
    close(
        step.stored_energy_change_j + step.viscous_heat_j,
        step.applied_work_j + step.background_work_j,
    );
    assert!(step.viscous_heat_j > 0.0);
}

#[test]
fn physical_regime_unsupported_drives_and_cancellation_still_refuse() {
    let source = sphere(20000.0);
    let body = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(1.0, 1.0);
    let make = || {
        CreepingSphere::new(&material, RADIUS, Gravity::ZERO, NativeVec3::ZERO, limits()).unwrap()
    };
    let state = initial(&body, NativeVec3::ZERO, NativeVec3::ZERO);
    assert!(matches!(
        make()
            .with_rotational_lift()
            .step_exact(&body, state, exact_controls(1.0), || false),
        Err(CreepingFlowError::Invalid(_))
    ));
    assert!(matches!(
        make().step_exact_with_wrench(
            &body,
            state,
            Wrench {
                force_world: NativeVec3::ZERO,
                torque_body: NativeVec3::new(0.0, 0.0, 1.0)
            },
            exact_controls(1.0),
            || false
        ),
        Err(CreepingFlowError::Invalid(_))
    ));
    let fast_drive = Wrench {
        force_world: NativeVec3::new(10.0 * 6.0 * core::f64::consts::PI * RADIUS, 0.0, 0.0),
        torque_body: NativeVec3::ZERO,
    };
    assert!(matches!(
        make().step_exact_with_wrench(&body, state, fast_drive, exact_controls(5.0), || false),
        Err(CreepingFlowError::Reynolds { .. })
    ));
    assert!(matches!(
        make().step_exact(&body, state, exact_controls(1.0), || true),
        Err(CreepingFlowError::Cancelled)
    ));
    for duration in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            make()
                .step_exact(&body, state, exact_controls(duration), || false)
                .is_err()
        );
    }
    let light_source = sphere(100.0);
    let light_body = light_source.principal_body(1e-12, || false).unwrap();
    let light_state = initial(&light_body, NativeVec3::ZERO, NativeVec3::ZERO);
    assert!(matches!(
        make().step_exact(&light_body, light_state, exact_controls(1.0), || false),
        Err(CreepingFlowError::Unsteady { .. })
    ));
}
