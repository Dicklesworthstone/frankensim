//! Original material and interface data drive exact native velocity relaxation.
#[path = "exact/flow.rs"]
mod flow;
use super::*;
use fs_mbd::{Pose, RigidBodyState, Vec3 as NativeVec3, Wrench};
use fs_mbd_material::creeping_flow::exact::ExactCreepingStepControls;
use fs_mbd_material::creeping_flow::{ExactCreepingAffineStep, ExactCreepingAffineStepControls};
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


fn affine_controls(
    duration_s: f64,
    minimum_force_variation_timescale_s: f64,
) -> ExactCreepingAffineStepControls {
    ExactCreepingAffineStepControls {
        step: exact_controls(duration_s),
        minimum_force_variation_timescale_s,
    }
}

fn affine_wrench(force_world: NativeVec3) -> Wrench {
    Wrench { force_world, torque_body: NativeVec3::ZERO }
}

fn affine_interpolate(first: NativeVec3, last: NativeVec3, fraction: f64) -> NativeVec3 {
    first.scale(1.0 - fraction).add(last.scale(fraction))
}

// Independent ordinary closed ODE solution at moderate relaxation exposure.
// The material body, fluid and Navier source still supply the actual step.
#[allow(clippy::too_many_arguments)]
fn affine_velocity(
    initial: NativeVec3,
    ambient: NativeVec3,
    first: NativeVec3,
    last: NativeVec3,
    net_weight: NativeVec3,
    mass: f64,
    resistance: f64,
    duration: f64,
    time: f64,
) -> NativeVec3 {
    let rate = resistance / mass;
    let slope = last.sub(first).scale(1.0 / (resistance * duration));
    let offset = first.add(net_weight).scale(1.0 / resistance)
        .sub(slope.scale(1.0 / rate));
    ambient.add(offset).add(slope.scale(time))
        .add(initial.sub(ambient).sub(offset).scale((-rate * time).exp()))
}

fn affine_power_integral(mut power: impl FnMut(f64) -> f64, duration: f64) -> f64 {
    let intervals = 8192;
    let dt = duration / f64::from(intervals);
    let mut sum = power(0.0) + power(duration);
    let mut correction = 0.0;
    for i in 1..intervals {
        let term = power(f64::from(i) * dt) * if i % 2 == 0 { 2.0 } else { 4.0 };
        let adjusted = term - correction;
        let next = sum + adjusted;
        correction = (next - sum) - adjusted;
        sum = next;
    }
    sum * dt / 3.0
}

fn close_integral(actual: f64, expected: f64) {
    assert!(expected.abs() > 1e-20, "fixture needs a nonzero power integral");
    assert!(
        (actual - expected).abs() <= 2e-11 * expected.abs(),
        "{actual} != {expected}",
    );
}

#[test]
fn material_and_navier_sources_drive_affine_motion_and_independent_work_heat() {
    let ambient = NativeVec3::new(0.03, -0.02, 0.01);
    let velocity = ambient.add(NativeVec3::new(0.03, -0.01, 0.02));
    let omega = NativeVec3::new(0.0, 0.0, 0.4);
    let first = NativeVec3::new(0.004, -0.002, 0.001);
    let last = NativeVec3::new(-0.003, 0.003, 0.002);
    let duration = 3.0;
    for (density, fluid_density, viscosity, slip) in [
        (20000.0, 1.0, 1.0, 0.0),
        (40000.0, 2.0, 2.0, 0.0),
        (20000.0, 1.0, 1.0, RADIUS),
    ] {
        let source = sphere(density);
        let body = source.principal_body(1e-12, || false).unwrap();
        let material = fluid(fluid_density, viscosity);
        let interface = super::navier_slip_source::constant(
            source.regions()[0].material().material(), slip, 400.0,
        );
        let reservoir = CreepingSphere::new(&material, RADIUS, gravity(), ambient, limits())
            .unwrap()
            .with_navier_slip(NavierSlipBoundary::new(&interface).unwrap())
            .unwrap();
        let state = initial(&body, velocity, omega);
        let result: ExactCreepingAffineStep = reservoir.step_exact_with_affine_force(
            &body, state, affine_wrench(first), last, affine_controls(duration, duration), || false,
        ).unwrap();
        let step = result.flow;
        let mass = body.properties().mass();
        let volume = 4.0 * core::f64::consts::PI * RADIUS.powi(3) / 3.0;
        let slip_ratio = slip / RADIUS;
        let c = 6.0 * core::f64::consts::PI * viscosity * RADIUS
            * (1.0 + 2.0 * slip_ratio) / (1.0 + 3.0 * slip_ratio);
        let c_rot = 8.0 * core::f64::consts::PI * viscosity * RADIUS.powi(3)
            / (1.0 + 3.0 * slip_ratio);
        let spin_rate = c_rot / body.properties().principal_inertia_body().x;
        let weight = gravity().acceleration_world().scale(mass - fluid_density * volume);
        let expected_v = |time| affine_velocity(
            velocity, ambient, first, last, weight, mass, c, duration, time,
        );
        close_vec(
            step.dynamics.state_after.center_of_mass_velocity_world(body.properties()).unwrap(),
            expected_v(duration),
        );
        close_vec(
            step.dynamics.displacement_world,
            NativeVec3::new(
                affine_power_integral(|t| expected_v(t).x, duration),
                affine_power_integral(|t| expected_v(t).y, duration),
                affine_power_integral(|t| expected_v(t).z, duration),
            ),
        );
        close_vec(
            body.properties().angular_velocity_body_checked(
                step.dynamics.state_after.angular_momentum_body(),
            ).unwrap(),
            omega.scale((-spin_rate * duration).exp()),
        );
        let external_work = affine_power_integral(
            |t| affine_interpolate(first, last, t / duration).dot(expected_v(t)), duration,
        );
        let independent_heat = affine_power_integral(
            |t| c * expected_v(t).sub(ambient).norm_squared()
                + c_rot * omega.norm_squared() * (-2.0 * spin_rate * t).exp(),
            duration,
        );
        close_integral(step.applied_work_j, external_work);
        close_integral(step.viscous_heat_j, independent_heat);
        let endpoint_mean_force_work = first.add(last).scale(0.5).dot(step.dynamics.displacement_world);
        assert!((external_work - endpoint_mean_force_work).abs() > 1e-3 * external_work.abs());
        close(step.stored_energy_change_j + step.viscous_heat_j,
            step.background_work_j + step.applied_work_j);
        assert!(step.energy_defect_j.abs() < 1e-11);
        assert!(result.maximum_translation_reynolds_bound <= limits().maximum_reynolds);
        close(result.drive_diffusion_ratio,
            RADIUS.powi(2) * fluid_density / (viscosity * duration));
        assert!(core::ptr::eq(reservoir.source(), &material));
        assert!(core::ptr::eq(reservoir.navier_slip().unwrap().source(), &interface));
    }
}

#[test]
fn affine_material_admission_refuses_an_interior_reynolds_peak_between_zero_slips() {
    let source = sphere(20000.0);
    let body = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(1.0, 1.0);
    let reservoir = CreepingSphere::new(
        &material, RADIUS, Gravity::ZERO, NativeVec3::ZERO, limits(),
    ).unwrap();
    let c = 6.0 * core::f64::consts::PI * RADIUS;
    let c_rot = 8.0 * core::f64::consts::PI * RADIUS.powi(3);
    let duration = body.properties().mass() / c;
    let velocity_scale = 20.0;
    let first = NativeVec3::new(c * velocity_scale, 0.0, 0.0);
    let last = first.scale(2.0 - core::f64::consts::E);
    let state = initial(&body, NativeVec3::ZERO, NativeVec3::ZERO);
    let native = fs_mbd::linear_drag::IsotropicLinearDrag::new(
        body.properties(), Gravity::ZERO, NativeVec3::ZERO, c, c_rot,
    ).unwrap().step_affine_force(state, affine_wrench(first), last, duration).unwrap();
    // Both actual endpoint loads admit, yet v(s)=V[e+(1-e)s-e*exp(-s)]
    // exceeds the same source's Reynolds ceiling in the interval.
    let before = reservoir.load(&body, state, || false).unwrap();
    let after = reservoir.load(&body, native.state_after, || false).unwrap();
    assert_eq!(before.translation_reynolds, 0.0);
    assert!(after.translation_reynolds < 1e-12);
    let peak_time_fraction = -(1.0 - (-1.0_f64).exp()).ln();
    let peak_speed = velocity_scale
        * (1.0 + (1.0 - core::f64::consts::E) * peak_time_fraction);
    let interior_reynolds = 2.0 * RADIUS * peak_speed;
    assert!(interior_reynolds > limits().maximum_reynolds);
    match reservoir.step_exact_with_affine_force(
        &body, state, affine_wrench(first), last, affine_controls(duration, duration), || false,
    ) {
        Err(CreepingFlowError::Reynolds { translation, rotation, maximum }) => {
            assert!(translation > interior_reynolds);
            close(translation, 2.0 * RADIUS * velocity_scale * -(-1.0_f64).exp_m1());
            assert_eq!(rotation, 0.0);
            assert_eq!(maximum, limits().maximum_reynolds);
        }
        result => panic!("interior-invalid affine drive must refuse: {result:?}"),
    }
}

#[test]
fn equal_force_retains_sharper_held_admission_and_the_identical_receipt() {
    let source = sphere(20000.0);
    let body = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(1.0, 1.0);
    let reservoir = CreepingSphere::new(
        &material, RADIUS, Gravity::ZERO, NativeVec3::ZERO, limits(),
    ).unwrap();
    let c = 6.0 * core::f64::consts::PI * RADIUS;
    let force = NativeVec3::new(6.0 * c, 0.0, 0.0);
    let state = initial(&body, NativeVec3::new(-2.0, 0.0, 0.0), NativeVec3::ZERO);
    let duration = -0.75_f64.ln() * body.properties().mass() / c;
    // The exact path runs from -2 to zero; the general triangle bound is 3.
    // Only the sharper original held admission fits the speed ceiling of 2.5.
    let held = reservoir.step_exact_with_wrench(
        &body, state, affine_wrench(force), exact_controls(duration), || false,
    ).unwrap();
    let affine = reservoir.step_exact_with_affine_force(
        &body, state, affine_wrench(force), force, affine_controls(duration, 0.0), || false,
    ).unwrap();
    assert_eq!(affine.flow, held);
    assert_eq!(affine.drive_diffusion_ratio, 0.0);
    close(affine.maximum_relative_speed_bound_m_s, 2.0);
    close(affine.maximum_translation_reynolds_bound, 0.04);
}

#[test]
fn source_affine_intervals_compose_and_late_cancellation_publishes_no_candidate() {
    let source = sphere(20000.0);
    let body = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(1.0, 1.0);
    let ambient = NativeVec3::new(0.03, -0.02, 0.01);
    let reservoir = CreepingSphere::new(&material, RADIUS, gravity(), ambient, limits()).unwrap();
    let state = initial(
        &body, ambient.add(NativeVec3::new(0.01, 0.03, -0.02)),
        NativeVec3::new(0.0, 0.0, 0.4),
    );
    let first = NativeVec3::new(0.004, -0.002, 0.001);
    let last = NativeVec3::new(-0.003, 0.003, 0.002);
    let duration = 3.0;
    let split = 0.7;
    let middle = affine_interpolate(first, last, split / duration);
    let mut polls = 0;
    let whole = reservoir.step_exact_with_affine_force(
        &body, state, affine_wrench(first), last, affine_controls(duration, duration),
        || { polls += 1; false },
    ).unwrap().flow;
    let a = reservoir.step_exact_with_affine_force(
        &body, state, affine_wrench(first), middle, affine_controls(split, duration), || false,
    ).unwrap().flow;
    let b = reservoir.step_exact_with_affine_force(
        &body, a.dynamics.state_after, affine_wrench(middle), last,
        affine_controls(duration - split, duration), || false,
    ).unwrap().flow;
    close_vec(b.dynamics.state_after.linear_momentum_world(),
        whole.dynamics.state_after.linear_momentum_world());
    close_vec(b.dynamics.state_after.pose().position_world(),
        whole.dynamics.state_after.pose().position_world());
    close_integral(a.viscous_heat_j + b.viscous_heat_j, whole.viscous_heat_j);
    close_integral(a.applied_work_j + b.applied_work_j, whole.applied_work_j);
    close_integral(a.background_work_j + b.background_work_j, whole.background_work_j);
    let mut visits = 0;
    assert!(polls > 1);
    assert!(matches!(reservoir.step_exact_with_affine_force(
        &body, state, affine_wrench(first), last, affine_controls(duration, duration),
        || { visits += 1; visits == polls },
    ), Err(CreepingFlowError::Cancelled)));
    assert_eq!(visits, polls);
    assert!(matches!(reservoir.with_rotational_lift().step_exact_with_affine_force(
        &body, state, affine_wrench(first), last, affine_controls(duration, duration), || false,
    ), Err(CreepingFlowError::Invalid(_))));
}


#[test]
fn physical_force_timescale_is_required_and_survives_numerical_subdivision() {
    let source = sphere(20000.0);
    let body = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(1.0, 1.0);
    let reservoir = CreepingSphere::new(
        &material, RADIUS, Gravity::ZERO, NativeVec3::ZERO, limits(),
    ).unwrap();
    let state = initial(&body, NativeVec3::ZERO, NativeVec3::ZERO);
    let first = NativeVec3::new(0.002, 0.0, 0.0);
    let last = NativeVec3::new(0.0, 0.003, 0.0);
    let source_timescale = 1.0;
    let whole = reservoir.step_exact_with_affine_force(
        &body, state, affine_wrench(first), last,
        affine_controls(1.0, source_timescale), || false,
    ).unwrap();
    let tiny_dt = 1e-6;
    let part = reservoir.step_exact_with_affine_force(
        &body, state, affine_wrench(first),
        affine_interpolate(first, last, tiny_dt),
        affine_controls(tiny_dt, source_timescale), || false,
    ).unwrap();
    assert_eq!(part.drive_diffusion_ratio, whole.drive_diffusion_ratio);
    close(whole.drive_diffusion_ratio, 1e-4);
    // A dt-based gate would incorrectly refuse this refinement.
    assert!(RADIUS.powi(2) / tiny_dt > limits().maximum_diffusion_ratio);
    for missing in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(reservoir.step_exact_with_affine_force(
            &body, state, affine_wrench(first), last,
            affine_controls(1.0, missing), || false,
        ), Err(CreepingFlowError::Invalid(_))));
    }
    match reservoir.step_exact_with_affine_force(
        &body, state, affine_wrench(first), last,
        affine_controls(1.0, 1e-5), || false,
    ) {
        Err(CreepingFlowError::UnsteadyDrive { diffusion_ratio, maximum }) => {
            close(diffusion_ratio, 10.0);
            assert_eq!(maximum, limits().maximum_diffusion_ratio);
        }
        result => panic!("fast physical forcing must refuse quasi-steady drag: {result:?}"),
    }
}
