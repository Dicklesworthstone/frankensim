//! Actual authored flow reaches source expansion, native motion and heat ownership.
use super::*;
use fs_mbd_material::thermal::heating::motion::creeping_flow::HeatedCreepingAffineFlow;

fn ambient(first: NativeVec3, last: NativeVec3, scale: f64) -> HeatedCreepingAffineFlow {
    HeatedCreepingAffineFlow {
        ambient_velocity_before_world: first,
        ambient_velocity_after_world: last,
        minimum_flow_variation_timescale_s: scale,
    }
}

fn close_vector(actual: NativeVec3, expected: NativeVec3) {
    close(actual.x, expected.x);
    close(actual.y, expected.y);
    close(actual.z, expected.z);
}

#[test]
fn expanding_navier_sphere_integrates_actual_flow_pressure_and_all_power_accounts() {
    use fs_mbd_material::creeping_flow::navier_slip::NavierSlipBoundary;
    let alpha = 0.02;
    let rho = 0.002;
    let mu = 0.1;
    let duration = 4.0 * DT;
    let first = NativeVec3::new(0.01, -0.005, 0.005);
    let last = NativeVec3::new(-0.005, 0.01, -0.01);
    let ambient_first = NativeVec3::new(0.04, -0.01, 0.02);
    let ambient_last = NativeVec3::new(-0.02, 0.03, -0.01);
    let acceleration = ambient_last.sub(ambient_first).scale(1.0 / duration);
    let material = card(CP, alpha, 500.0);
    let source = reference_body(&material, VolumeMoments::sphere(RADIUS).unwrap());
    let chart = caloric(&material, 1e5, SolidLiquidPhase::Solid);
    let original = body(
        &source,
        &material,
        &chart,
        RATE,
        controls(),
        gravity(),
        twist(),
    );
    let reservoir = fluid(rho, mu);
    let interface = navier_slip_source::constant(material.id(), RADIUS, 400.0);
    // If a source window accidentally reuses the constructor bath, admission
    // and the resulting trajectory disagree with this authored history.
    let sphere = CreepingSphere::new(
        &reservoir,
        RADIUS,
        gravity(),
        NativeVec3::new(30.0, 10.0, 20.0),
        limits(),
    )
    .unwrap()
    .with_navier_slip(NavierSlipBoundary::new(&interface).unwrap())
    .unwrap();
    let (next, report) = original
        .advance_creeping_sphere_exact_affine_flow(
            &sphere,
            drive(first, last, 2.0 * duration),
            ambient(ambient_first, ambient_last, 3.0 * duration),
            exact_controls(),
            duration,
            config(),
            || Ok(()),
        )
        .unwrap();

    let mass = source.inertia().mass();
    let volume = |r: f64| 4.0 * std::f64::consts::PI * r.powi(3) / 3.0;
    let inertia = |r: f64| 0.4 * mass * r * r;
    let resistance = |r: f64| {
        6.0 * std::f64::consts::PI * mu * r * (1.0 + 2.0 * RADIUS / r) / (1.0 + 3.0 * RADIUS / r)
    };
    let mut velocity = NativeVec3::new(0.01, -0.02, 0.03);
    let mut position = NativeVec3::new(0.0, 0.0, -0.2);
    let mut spin = 0.2 * inertia(RADIUS);
    let mut applied_work = 0.0;
    let mut pressure_work = 0.0;
    let mut heat = 0.0;
    let mut background_work = 0.0;
    let mut parameter_work = 0.0;
    let mut maximum_reynolds: f64 = 0.0;
    for index in 0..4 {
        let time = f64::from(index) * DT;
        let radius = RADIUS * (1.0 + alpha * RATE * time);
        let next_radius = RADIUS * (1.0 + alpha * RATE * (time + DT));
        let c = resistance(radius);
        let rotary =
            8.0 * std::f64::consts::PI * mu * radius.powi(3) / (1.0 + 3.0 * RADIUS / radius);
        let rate = c / mass;
        let before = interpolate_force(first, last, time / duration);
        let after = interpolate_force(first, last, (time + DT) / duration);
        let u0 = interpolate_force(ambient_first, ambient_last, time / duration);
        let pressure = acceleration.scale(rho * volume(radius));
        let weight = gravity()
            .acceleration_world()
            .scale(mass - rho * volume(radius));
        // Independent absolute ODE: c*U(t), affine applied force and the held
        // displaced-fluid pressure are the inhomogeneous polynomial drive.
        let slope = after.sub(before).scale(1.0 / (c * DT)).add(acceleration);
        let offset = before
            .add(weight)
            .add(pressure)
            .scale(1.0 / c)
            .sub(slope.scale(1.0 / rate));
        let amplitude = velocity.sub(u0).sub(offset);
        let path = |t: f64| {
            u0.add(offset)
                .add(slope.scale(t))
                .add(amplitude.scale((-rate * t).exp()))
        };
        let bath = |t: f64| u0.add(acceleration.scale(t));
        applied_work += integrate_power(
            |t| interpolate_force(before, after, t / DT).dot(path(t)),
            DT,
        );
        let interval_pressure_work = integrate_power(|t| pressure.dot(path(t)), DT);
        pressure_work += interval_pressure_work;
        heat += integrate_power(|t| c * path(t).sub(bath(t)).norm_squared(), DT)
            + 0.5 * spin * spin / inertia(radius)
                * -(-2.0 * DT * rotary / inertia(radius)).exp_m1();
        background_work += interval_pressure_work
            + integrate_power(|t| -c * path(t).sub(bath(t)).dot(bath(t)), DT);
        let drive_bound = speed(
            before
                .add(weight)
                .add(pressure)
                .scale(1.0 / mass)
                .sub(acceleration),
        )
        .max(speed(
            after
                .add(weight)
                .add(pressure)
                .scale(1.0 / mass)
                .sub(acceleration),
        ));
        let initial_slip = speed(velocity.sub(u0));
        let bound = initial_slip
            .max((-rate * DT).exp() * initial_slip + drive_bound * -(-rate * DT).exp_m1() / rate);
        maximum_reynolds = maximum_reynolds.max(2.0 * radius * rho * bound / mu);
        position = position
            .add(u0.add(offset).scale(DT))
            .add(slope.scale(0.5 * DT * DT))
            .add(amplitude.scale(-(-rate * DT).exp_m1() / rate));
        velocity = path(DT);
        spin *= (-DT * rotary / inertia(radius)).exp();
        parameter_work += 0.5 * spin * spin * (1.0 / inertia(next_radius) - 1.0 / inertia(radius))
            + rho
                * (volume(next_radius) - volume(radius))
                * gravity().acceleration_world().dot(position);
        maximum_reynolds =
            maximum_reynolds.max(2.0 * next_radius * rho * speed(velocity.sub(bath(DT))) / mu);
    }
    close_vector(
        next.state()
            .center_of_mass_velocity_world(next.principal_frame().properties())
            .unwrap(),
        velocity,
    );
    close_vector(next.state().pose().position_world(), position);
    close(
        next.state()
            .pose()
            .orientation()
            .rotate_body_to_world(next.state().angular_momentum_body())
            .z,
        spin,
    );
    assert!(applied_work.abs() > 1e-9 && pressure_work.abs() > 1e-9 && heat > 0.0);
    close(report.applied_work_j / applied_work, 1.0);
    close(report.ambient_pressure_work_j / pressure_work, 1.0);
    close(report.viscous_heat_j / heat, 1.0);
    close(report.background_work_j / background_work, 1.0);
    close(report.parameter_work_j, parameter_work);
    close(report.maximum_translation_reynolds_bound, maximum_reynolds);
    let final_radius = RADIUS * (1.0 + alpha * RATE * duration);
    close(
        report.maximum_drive_diffusion_ratio,
        final_radius.powi(2) * rho / (mu * 2.0 * duration),
    );
    close(
        report.maximum_flow_diffusion_ratio,
        final_radius.powi(2) * rho / (mu * 3.0 * duration),
    );
    close_vector(
        report.load_after.drag_force_world_n,
        velocity.sub(ambient_last).scale(-resistance(final_radius)),
    );
    assert_eq!(report.ambient_velocity_before_world, ambient_first);
    assert_eq!(report.ambient_velocity_after_world, ambient_last);
    close(report.thermal.external_heat_j, mass * CP * RATE * duration);
    close(report.combined_residual_j, 0.0);
    assert_eq!(next.configuration().inertia().mass(), mass);
    assert_eq!(sphere.navier_slip().unwrap().slip_length_m(), RADIUS);
    assert!(core::ptr::eq(sphere.source(), &reservoir));
    assert!(core::ptr::eq(
        sphere.navier_slip().unwrap().source(),
        &interface
    ));
}

#[test]
fn flow_history_composes_on_actual_epochs_and_held_flow_retains_original_parity() {
    let material = card(CP, 0.0, 500.0);
    let source = reference_body(&material, VolumeMoments::sphere(RADIUS).unwrap());
    let chart = caloric(&material, 1e5, SolidLiquidPhase::Solid);
    let mut coarse_controls = controls();
    coarse_controls.maximum_step_s = 2.0 * DT;
    let original = body(
        &source,
        &material,
        &chart,
        RATE,
        controls(),
        gravity(),
        twist(),
    );
    let coarse = body(
        &source,
        &material,
        &chart,
        RATE,
        coarse_controls,
        gravity(),
        twist(),
    );
    let reservoir = fluid(0.002, 0.1);
    let sphere = CreepingSphere::new(&reservoir, RADIUS, gravity(), flow(), limits()).unwrap();
    let duration = 4.0 * DT;
    let first = NativeVec3::new(0.002, -0.001, 0.001);
    let last = NativeVec3::new(-0.001, 0.002, -0.002);
    let u0 = NativeVec3::new(0.04, -0.01, 0.02);
    let u1 = NativeVec3::new(-0.02, 0.03, -0.01);
    let force = drive(first, last, duration);
    let flow_history = ambient(u0, u1, duration);
    let (whole, audit) = original
        .advance_creeping_sphere_exact_affine_flow(
            &sphere,
            force,
            flow_history,
            exact_controls(),
            duration,
            config(),
            || Ok(()),
        )
        .unwrap();
    let (partitioned, partition_audit) = coarse
        .advance_creeping_sphere_exact_affine_flow(
            &sphere,
            force,
            flow_history,
            exact_controls(),
            duration,
            config(),
            || Ok(()),
        )
        .unwrap();
    let middle_force = interpolate_force(first, last, 0.5);
    let middle_flow = interpolate_force(u0, u1, 0.5);
    let (middle, first_audit) = original
        .advance_creeping_sphere_exact_affine_flow(
            &sphere,
            drive(first, middle_force, duration),
            ambient(u0, middle_flow, duration),
            exact_controls(),
            2.0 * DT,
            config(),
            || Ok(()),
        )
        .unwrap();
    let (continued, last_audit) = middle
        .advance_creeping_sphere_exact_affine_flow(
            &sphere,
            drive(middle_force, last, duration),
            ambient(first_audit.ambient_velocity_after_world, u1, duration),
            exact_controls(),
            duration,
            config(),
            || Ok(()),
        )
        .unwrap();
    for candidate in [&partitioned, &continued] {
        close_vector(
            candidate.state().linear_momentum_world(),
            whole.state().linear_momentum_world(),
        );
        close_vector(
            candidate.state().pose().position_world(),
            whole.state().pose().position_world(),
        );
    }
    for (total, coarse_value, pieces) in [
        (
            audit.applied_work_j,
            partition_audit.applied_work_j,
            first_audit.applied_work_j + last_audit.applied_work_j,
        ),
        (
            audit.background_work_j,
            partition_audit.background_work_j,
            first_audit.background_work_j + last_audit.background_work_j,
        ),
        (
            audit.ambient_pressure_work_j,
            partition_audit.ambient_pressure_work_j,
            first_audit.ambient_pressure_work_j + last_audit.ambient_pressure_work_j,
        ),
        (
            audit.viscous_heat_j,
            partition_audit.viscous_heat_j,
            first_audit.viscous_heat_j + last_audit.viscous_heat_j,
        ),
    ] {
        close(coarse_value / total, 1.0);
        close(pieces / total, 1.0);
    }
    assert_eq!(last_audit.ambient_velocity_before_world, middle_flow);
    assert_eq!(last_audit.ambient_velocity_after_world, u1);
    assert_eq!(
        partition_audit.maximum_flow_diffusion_ratio,
        audit.maximum_flow_diffusion_ratio
    );
    assert_eq!(
        last_audit.maximum_flow_diffusion_ratio,
        audit.maximum_flow_diffusion_ratio
    );

    let (held, held_audit) = original
        .advance_creeping_sphere_exact_affine(
            &sphere,
            force,
            exact_controls(),
            duration,
            config(),
            || Ok(()),
        )
        .unwrap();
    let held_end = NativeVec3::new(flow().x, -0.0, flow().z);
    let (compatible, compatible_audit) = original
        .advance_creeping_sphere_exact_affine_flow(
            &sphere,
            force,
            ambient(flow(), held_end, f64::NAN),
            exact_controls(),
            duration,
            config(),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(compatible.state(), held.state());
    assert_eq!(compatible_audit, held_audit);
    assert_eq!(compatible_audit.ambient_pressure_work_j, 0.0);
    assert_eq!(compatible_audit.maximum_flow_diffusion_ratio, 0.0);
    assert_eq!(
        compatible_audit.ambient_velocity_after_world.y.to_bits(),
        (-0.0_f64).to_bits()
    );
}

#[test]
fn changing_flow_and_generated_heat_publish_together_or_roll_back_together() {
    let material = card(CP, 0.02, 500.0);
    let source = reference_body(&material, VolumeMoments::sphere(RADIUS).unwrap());
    let chart = caloric(&material, 1e5, SolidLiquidPhase::Solid);
    let direct = body(
        &source,
        &material,
        &chart,
        RATE,
        controls(),
        gravity(),
        twist(),
    );
    let coupled = body(
        &source,
        &material,
        &chart,
        0.0,
        controls(),
        gravity(),
        twist(),
    );
    let reservoir = fluid(0.002, 0.1);
    let sphere = CreepingSphere::new(&reservoir, RADIUS, gravity(), flow(), limits()).unwrap();
    let force = drive(NativeVec3::ZERO, NativeVec3::ZERO, 0.0);
    let flow_history = ambient(flow(), NativeVec3::new(-0.01, 0.03, -0.02), 4.0 * DT);
    let input = PrescribedHeat {
        watts: source.inertia().mass() * CP * RATE,
        time_s: 0.0,
        supplied_j: 0.0,
        refuse_after_s: None,
    };
    let mut polls = 0;
    let (expected, expected_audit) = direct
        .advance_creeping_sphere_exact_affine_flow(
            &sphere,
            force,
            flow_history,
            exact_controls(),
            4.0 * DT,
            config(),
            || {
                polls += 1;
                Ok(())
            },
        )
        .unwrap();
    let (actual, generated, audit) = coupled
        .advance_creeping_sphere_exact_affine_flow_with(
            &input,
            &sphere,
            force,
            flow_history,
            exact_controls(),
            4.0 * DT,
            config(),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(actual.state(), expected.state());
    assert_eq!(actual.configuration(), expected.configuration());
    close(generated.supplied_j, audit.thermal.external_heat_j);
    close(
        audit.ambient_pressure_work_j,
        expected_audit.ambient_pressure_work_j,
    );
    close(audit.viscous_heat_j, expected_audit.viscous_heat_j);
    close(audit.combined_residual_j, 0.0);
    let failure = PrescribedHeat {
        refuse_after_s: Some(DT),
        ..input.clone()
    };
    assert!(
        coupled
            .advance_creeping_sphere_exact_affine_flow_with(
                &failure,
                &sphere,
                force,
                flow_history,
                exact_controls(),
                4.0 * DT,
                config(),
                || Ok(()),
            )
            .is_err()
    );
    assert_eq!(failure.time_s, 0.0);
    assert_eq!(failure.supplied_j, 0.0);
    let mut visits = 0;
    assert!(matches!(
        direct.advance_creeping_sphere_exact_affine_flow(
            &sphere,
            force,
            flow_history,
            exact_controls(),
            4.0 * DT,
            config(),
            || {
                visits += 1;
                if visits == polls {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            },
        ),
        Err(HeatedMotionError::Heating(
            fs_mbd_material::thermal::heating::HeatingError::Thermal(DepositionError::Cancelled)
        ))
    ));
    assert_eq!(visits, polls);
    assert_eq!(coupled.thermal().time_s(), 0.0);
    assert_eq!(direct.thermal().time_s(), 0.0);
    assert_eq!(coupled.state(), direct.state());
}

#[test]
fn expanded_endpoints_recheck_actual_flow_slip_and_original_physical_timescale() {
    let material = card(CP, 0.02, 500.0);
    let source = reference_body(&material, VolumeMoments::sphere(RADIUS).unwrap());
    let chart = caloric(&material, 1e5, SolidLiquidPhase::Solid);
    let original = body(
        &source,
        &material,
        &chart,
        RATE,
        controls(),
        gravity(),
        twist(),
    );
    let rho = 0.002;
    let mu = 0.001;
    let reservoir = fluid(rho, mu);
    let sphere = CreepingSphere::new(&reservoir, RADIUS, gravity(), flow(), limits()).unwrap();
    let source_scale = RADIUS.powi(2) * rho / (mu * 0.048);
    let first = flow();
    let last = flow().add(NativeVec3::new(1e-6, 0.0, 0.0));
    let force = drive(NativeVec3::ZERO, NativeVec3::ZERO, 0.0);
    let (_, prefix) = original
        .advance_creeping_sphere_exact_affine_flow(
            &sphere,
            force,
            ambient(first, interpolate_force(first, last, 0.5), source_scale),
            exact_controls(),
            2.0 * DT,
            config(),
            || Ok(()),
        )
        .unwrap();
    close(
        prefix.maximum_flow_diffusion_ratio,
        0.048 * 1.02_f64.powi(2),
    );
    assert!(prefix.maximum_diffusion_ratio < limits().maximum_diffusion_ratio);
    let mut refined_controls = controls();
    refined_controls.maximum_step_s = DT / 2.0;
    let refined = body(
        &source,
        &material,
        &chart,
        RATE,
        refined_controls,
        gravity(),
        twist(),
    );
    for candidate in [&original, &refined] {
        assert!(matches!(
            candidate.advance_creeping_sphere_exact_affine_flow(
                &sphere,
                force,
                ambient(first, last, source_scale),
                exact_controls(),
                4.0 * DT,
                config(),
                || Ok(()),
            ),
            Err(HeatedMotionError::CreepingFlow(
                CreepingFlowError::UnsteadyDrive { .. }
            ))
        ));
        assert_eq!(candidate.thermal().time_s(), 0.0);
    }
    for scale in [0.0, f64::NAN] {
        assert!(matches!(
            original.advance_creeping_sphere_exact_affine_flow(
                &sphere,
                force,
                ambient(first, last, scale),
                exact_controls(),
                DT,
                config(),
                || Ok(()),
            ),
            Err(HeatedMotionError::Invalid(_))
        ));
    }

    // This exact interval passes at the starting radius. Only the updated
    // radius AND the actual accepted ambient endpoint expose the invalid slip.
    let initial_speed = 0.12458;
    let after_flow = NativeVec3::new(-0.0001, 0.0, 0.0);
    let resting_spin = Twist::new(Vec3::default(), Vec3::new(initial_speed, 0.0, 0.0));
    let endpoint_body = body(
        &source,
        &material,
        &chart,
        RATE,
        controls(),
        Gravity::ZERO,
        resting_spin,
    );
    let endpoint_sphere = CreepingSphere::new(
        &reservoir,
        RADIUS,
        Gravity::ZERO,
        NativeVec3::ZERO,
        limits(),
    )
    .unwrap();
    let rate = 6.0 * std::f64::consts::PI * mu * RADIUS / source.inertia().mass();
    let relative_drive = -(1.0 - rho / DENSITY) * after_flow.x / DT;
    let final_slip =
        initial_speed * (-rate * DT).exp() + relative_drive * -(-rate * DT).exp_m1() / rate;
    let re_factor = 2.0 * RADIUS * 1.01 * rho / mu;
    assert!(2.0 * RADIUS * rho * initial_speed / mu < limits().maximum_reynolds);
    assert!(re_factor * (final_slip + after_flow.x) < limits().maximum_reynolds);
    assert!(re_factor * final_slip > limits().maximum_reynolds);
    assert!(matches!(
        endpoint_body.advance_creeping_sphere_exact_affine_flow(
            &endpoint_sphere,
            force,
            ambient(NativeVec3::ZERO, after_flow, 1.0),
            exact_controls(),
            DT,
            config(),
            || Ok(()),
        ),
        Err(HeatedMotionError::CreepingFlow(
            CreepingFlowError::Reynolds { .. }
        ))
    ));
    assert_eq!(endpoint_body.thermal().time_s(), 0.0);
}
