//! Actual source/geometry and ambient history reach the native exact trajectory.
use super::*;
use fs_mbd_material::creeping_flow::ExactCreepingFlowStepControls;

fn flow_controls(duration: f64, source_timescale: f64) -> ExactCreepingFlowStepControls {
    ExactCreepingFlowStepControls {
        step: exact_controls(duration),
        minimum_force_variation_timescale_s: source_timescale,
        minimum_flow_variation_timescale_s: source_timescale,
    }
}

#[test]
fn sourced_affine_flow_retains_pressure_and_independent_force_bath_and_heat_integrals() {
    let ambient = NativeVec3::new(0.03, -0.02, 0.01);
    let ambient_after = NativeVec3::new(-0.01, 0.04, -0.02);
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
        let interface = super::super::navier_slip_source::constant(
            source.regions()[0].material().material(),
            slip,
            400.0,
        );
        let reservoir = CreepingSphere::new(&material, RADIUS, gravity(), ambient, limits())
            .unwrap()
            .with_navier_slip(NavierSlipBoundary::new(&interface).unwrap())
            .unwrap();
        let result = reservoir
            .step_exact_with_affine_force_and_flow(
                &body,
                initial(&body, velocity, omega),
                affine_wrench(first),
                last,
                ambient,
                ambient_after,
                flow_controls(duration, duration),
                || false,
            )
            .unwrap();
        let mass = body.properties().mass();
        let fluid_mass = fluid_density * 4.0 * core::f64::consts::PI * RADIUS.powi(3) / 3.0;
        let ratio = slip / RADIUS;
        let c = 6.0 * core::f64::consts::PI * viscosity * RADIUS * (1.0 + 2.0 * ratio)
            / (1.0 + 3.0 * ratio);
        let c_rot = 8.0 * core::f64::consts::PI * viscosity * RADIUS.powi(3) / (1.0 + 3.0 * ratio);
        let pressure = ambient_after.sub(ambient).scale(fluid_mass / duration);
        let spin_rate = c_rot / body.properties().principal_inertia_body().x;
        // Independent absolute ODE: c*(U(t)-U0) is an affine drive,
        // while displaced-fluid pressure is the separate held physical force.
        let v = |t| {
            affine_velocity(
                velocity,
                ambient,
                first.add(pressure),
                last.add(pressure).add(ambient_after.sub(ambient).scale(c)),
                gravity().acceleration_world().scale(mass - fluid_mass),
                mass,
                c,
                duration,
                t,
            )
        };
        let u = |t| affine_interpolate(ambient, ambient_after, t / duration);
        let step = result.flow;
        close_vec(
            step.dynamics
                .state_after
                .center_of_mass_velocity_world(body.properties())
                .unwrap(),
            v(duration),
        );
        close_vec(
            step.dynamics.displacement_world,
            NativeVec3::new(
                affine_power_integral(|t| v(t).x, duration),
                affine_power_integral(|t| v(t).y, duration),
                affine_power_integral(|t| v(t).z, duration),
            ),
        );
        close_vec(result.ambient_pressure_force_world_n, pressure);
        close_integral(
            result.ambient_pressure_work_j,
            affine_power_integral(|t| pressure.dot(v(t)), duration),
        );
        close_integral(
            step.applied_work_j,
            affine_power_integral(
                |t| affine_interpolate(first, last, t / duration).dot(v(t)),
                duration,
            ),
        );
        close_integral(
            step.background_work_j,
            affine_power_integral(
                |t| pressure.dot(v(t)) - c * v(t).sub(u(t)).dot(u(t)),
                duration,
            ),
        );
        close_integral(
            step.viscous_heat_j,
            affine_power_integral(
                |t| {
                    c * v(t).sub(u(t)).norm_squared()
                        + c_rot * omega.norm_squared() * (-2.0 * spin_rate * t).exp()
                },
                duration,
            ),
        );
        close_vec(
            step.load_after.drag_force_world_n,
            v(duration).sub(ambient_after).scale(-c),
        );
        close(
            step.stored_energy_change_j + step.viscous_heat_j,
            step.background_work_j + step.applied_work_j,
        );
        assert!(step.energy_defect_j.abs() < 1e-11);
        assert_eq!(result.ambient_velocity_after_world, ambient_after);
        close(
            result.flow_diffusion_ratio,
            RADIUS.powi(2) * fluid_density / (viscosity * duration),
        );
        assert_eq!(result.flow_diffusion_ratio, result.drive_diffusion_ratio);
        assert!(core::ptr::eq(reservoir.source(), &material));
        assert!(core::ptr::eq(
            reservoir.navier_slip().unwrap().source(),
            &interface
        ));
    }
}

#[test]
fn explicit_accepted_ambient_endpoints_compose_without_restarting_the_sphere_bath() {
    let source = sphere(20000.0);
    let body = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(1.0, 1.0);
    // The immutable constructor's bath deliberately differs from this history.
    let reservoir = CreepingSphere::new(
        &material,
        RADIUS,
        gravity(),
        NativeVec3::new(9.0, 8.0, 7.0),
        limits(),
    )
    .unwrap();
    let ambient = NativeVec3::new(0.03, -0.02, 0.01);
    let end = NativeVec3::new(-0.01, 0.04, -0.02);
    let state = initial(
        &body,
        ambient.add(NativeVec3::new(0.01, 0.02, -0.01)),
        NativeVec3::ZERO,
    );
    let first = NativeVec3::new(0.004, -0.002, 0.001);
    let last = NativeVec3::new(-0.003, 0.003, 0.002);
    let duration = 3.0;
    let split = 0.7;
    let middle = affine_interpolate(first, last, split / duration);
    let middle_flow = affine_interpolate(ambient, end, split / duration);
    let mut polls = 0;
    let whole = reservoir
        .step_exact_with_affine_force_and_flow(
            &body,
            state,
            affine_wrench(first),
            last,
            ambient,
            end,
            flow_controls(duration, duration),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let a = reservoir
        .step_exact_with_affine_force_and_flow(
            &body,
            state,
            affine_wrench(first),
            middle,
            ambient,
            middle_flow,
            flow_controls(split, duration),
            || false,
        )
        .unwrap();
    let b = reservoir
        .step_exact_with_affine_force_and_flow(
            &body,
            a.flow.dynamics.state_after,
            affine_wrench(middle),
            last,
            a.ambient_velocity_after_world,
            end,
            flow_controls(duration - split, duration),
            || false,
        )
        .unwrap();
    close_vec(
        b.flow.dynamics.state_after.linear_momentum_world(),
        whole.flow.dynamics.state_after.linear_momentum_world(),
    );
    close_vec(
        b.flow.dynamics.state_after.pose().position_world(),
        whole.flow.dynamics.state_after.pose().position_world(),
    );
    close_integral(
        a.flow.viscous_heat_j + b.flow.viscous_heat_j,
        whole.flow.viscous_heat_j,
    );
    close_integral(
        a.flow.applied_work_j + b.flow.applied_work_j,
        whole.flow.applied_work_j,
    );
    close_integral(
        a.flow.background_work_j + b.flow.background_work_j,
        whole.flow.background_work_j,
    );
    close_integral(
        a.ambient_pressure_work_j + b.ambient_pressure_work_j,
        whole.ambient_pressure_work_j,
    );
    assert_eq!(a.flow_diffusion_ratio, whole.flow_diffusion_ratio);
    assert_eq!(b.flow_diffusion_ratio, whole.flow_diffusion_ratio);
    let mut visits = 0;
    assert!(matches!(
        reservoir.step_exact_with_affine_force_and_flow(
            &body,
            state,
            affine_wrench(first),
            last,
            ambient,
            end,
            flow_controls(duration, duration),
            || {
                visits += 1;
                visits == polls
            },
        ),
        Err(CreepingFlowError::Cancelled)
    ));
    assert_eq!(visits, polls);
}

#[test]
fn source_flow_timescale_survives_refinement_and_held_flow_keeps_original_parity() {
    let source = sphere(20000.0);
    let body = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(1.0, 1.0);
    let reservoir =
        CreepingSphere::new(&material, RADIUS, Gravity::ZERO, NativeVec3::ZERO, limits()).unwrap();
    let state = initial(&body, NativeVec3::ZERO, NativeVec3::ZERO);
    let first = NativeVec3::new(0.002, 0.0, 0.0);
    let last = NativeVec3::new(0.0, 0.003, 0.0);
    let end = NativeVec3::new(0.01, -0.02, 0.0);
    let whole = reservoir
        .step_exact_with_affine_force_and_flow(
            &body,
            state,
            affine_wrench(first),
            last,
            NativeVec3::ZERO,
            end,
            flow_controls(1.0, 1.0),
            || false,
        )
        .unwrap();
    let dt = 1e-6;
    let part = reservoir
        .step_exact_with_affine_force_and_flow(
            &body,
            state,
            affine_wrench(first),
            affine_interpolate(first, last, dt),
            NativeVec3::ZERO,
            end.scale(dt),
            flow_controls(dt, 1.0),
            || false,
        )
        .unwrap();
    assert_eq!(part.flow_diffusion_ratio, whole.flow_diffusion_ratio);
    close(part.flow_diffusion_ratio, 1e-4);
    for scale in [0.0, f64::NAN] {
        let mut controls = flow_controls(1.0, 1.0);
        controls.minimum_flow_variation_timescale_s = scale;
        assert!(matches!(
            reservoir.step_exact_with_affine_force_and_flow(
                &body,
                state,
                affine_wrench(first),
                last,
                NativeVec3::ZERO,
                end,
                controls,
                || false,
            ),
            Err(CreepingFlowError::Invalid(_))
        ));
    }
    let mut fast = flow_controls(dt, 1.0);
    fast.minimum_flow_variation_timescale_s = 1e-5;
    assert!(matches!(
        reservoir.step_exact_with_affine_force_and_flow(
            &body,
            state,
            affine_wrench(first),
            affine_interpolate(first, last, dt),
            NativeVec3::ZERO,
            end.scale(dt),
            fast,
            || false,
        ),
        Err(CreepingFlowError::UnsteadyDrive { .. })
    ));
    let held = reservoir
        .step_exact_with_affine_force(
            &body,
            state,
            affine_wrench(first),
            last,
            affine_controls(1.0, 1.0),
            || false,
        )
        .unwrap();
    let mut controls = flow_controls(1.0, 1.0);
    controls.minimum_flow_variation_timescale_s = f64::NAN;
    let compatible = reservoir
        .step_exact_with_affine_force_and_flow(
            &body,
            state,
            affine_wrench(first),
            last,
            NativeVec3::ZERO,
            NativeVec3::ZERO,
            controls,
            || false,
        )
        .unwrap();
    assert_eq!(compatible.flow, held.flow);
    assert_eq!(compatible.ambient_pressure_work_j, 0.0);
    let light_source = sphere(1.0);
    let light = light_source.principal_body(1e-12, || false).unwrap();
    assert!(matches!(
        reservoir.step_exact_with_affine_force_and_flow(
            &light,
            initial(&light, NativeVec3::ZERO, NativeVec3::ZERO),
            Wrench::ZERO,
            NativeVec3::ZERO,
            NativeVec3::ZERO,
            end,
            flow_controls(1.0, 1.0),
            || false,
        ),
        Err(CreepingFlowError::Unsteady { .. })
    ));
}

#[test]
fn complete_moving_flow_admission_catches_an_interior_slip_peak() {
    let source = sphere(20000.0);
    let body = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(1.0, 1.0);
    let reservoir =
        CreepingSphere::new(&material, RADIUS, Gravity::ZERO, NativeVec3::ZERO, limits()).unwrap();
    let c = 6.0 * core::f64::consts::PI * RADIUS;
    let duration = body.properties().mass() / c;
    let end = NativeVec3::new(1.0, 0.0, 0.0);
    let fluid_mass = 4.0 * core::f64::consts::PI * RADIUS.powi(3) / 3.0;
    let pressure = end.scale(fluid_mass / duration);
    let acceleration_force = end.scale(body.properties().mass() / duration);
    let drive = NativeVec3::new(20.0 * c, 0.0, 0.0);
    let first = drive.add(acceleration_force).sub(pressure);
    let last = drive
        .scale(2.0 - core::f64::consts::E)
        .add(acceleration_force)
        .sub(pressure);
    // The relative solution is V[e+(1-e)s-e exp(-s)], zero at both ends.
    let peak_fraction = -(1.0 - (-1.0_f64).exp()).ln();
    let peak_reynolds = 2.0 * RADIUS * 20.0 * (1.0 + (1.0 - core::f64::consts::E) * peak_fraction);
    assert!(peak_reynolds > limits().maximum_reynolds);
    match reservoir.step_exact_with_affine_force_and_flow(
        &body,
        initial(&body, NativeVec3::ZERO, NativeVec3::ZERO),
        affine_wrench(first),
        last,
        NativeVec3::ZERO,
        end,
        flow_controls(duration, duration),
        || false,
    ) {
        Err(CreepingFlowError::Reynolds { translation, .. }) => {
            assert!(translation > peak_reynolds);
            close(translation, 2.0 * RADIUS * 20.0 * -(-1.0_f64).exp_m1());
        }
        result => panic!("interior-invalid changing flow must refuse: {result:?}"),
    }
}
