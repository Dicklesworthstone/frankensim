//! Analytic frozen-radius motion on the existing accepted thermal timeline.
use super::*;
use fs_mbd_material::thermal::heating::motion::creeping_flow::{
    HeatedCreepingAffineForce, HeatedExactCreepingControls,
};
use fs_mbd::Vec3 as NativeVec3;

#[path = "exact/flow.rs"]
mod flow;

fn exact_controls() -> HeatedExactCreepingControls {
    HeatedExactCreepingControls {
        maximum_step_energy_defect_j: 1e-13,
    }
}

#[test]
fn stiff_exact_motion_advances_without_heating_the_finite_solid_with_fluid_loss() {
    let material = card(CP, 0.0, 500.0);
    let source = reference_body(&material, VolumeMoments::sphere(RADIUS).unwrap());
    let chart = caloric(&material, 1e5, SolidLiquidPhase::Solid);
    let original = body(
        &source,
        &material,
        &chart,
        0.0,
        controls(),
        gravity(),
        twist(),
    );
    let rho = 0.002;
    let mu = 1.0;
    let reservoir = fluid(rho, mu);
    let sphere = CreepingSphere::new(&reservoir, RADIUS, gravity(), flow(), limits()).unwrap();
    assert!(
        original
            .advance_creeping_sphere(&sphere, flow_controls(), DT, config(), || Ok(()))
            .is_err()
    );
    let (next, report) = original
        .advance_creeping_sphere_exact(&sphere, exact_controls(), DT, config(), || Ok(()))
        .unwrap();
    let mass = source.inertia().mass();
    let inertia = 0.4 * mass * RADIUS.powi(2);
    let c = 6.0 * std::f64::consts::PI * mu * RADIUS;
    let c_rot = 8.0 * std::f64::consts::PI * mu * RADIUS.powi(3);
    let tau = mass / c;
    let tau_rot = inertia / c_rot;
    let fluid_mass = rho * 4.0 * std::f64::consts::PI * RADIUS.powi(3) / 3.0;
    let terminal = flow().add(
        gravity()
            .acceleration_world()
            .scale((mass - fluid_mass) / c),
    );
    let initial = fs_mbd::Vec3::new(0.01, -0.02, 0.03);
    let deviation = initial.sub(terminal);
    let expected_v = terminal.add(deviation.scale((-DT / tau).exp()));
    let expected_dx = terminal
        .scale(DT)
        .add(deviation.scale(tau * -(-DT / tau).exp_m1()));
    let velocity = next
        .state()
        .center_of_mass_velocity_world(next.principal_frame().properties())
        .unwrap();
    let position = next.state().pose().position_world();
    for (actual, expected) in [
        (velocity.x, expected_v.x),
        (velocity.y, expected_v.y),
        (velocity.z, expected_v.z),
        (position.x, expected_dx.x),
        (position.y, expected_dx.y),
        (position.z, -0.2 + expected_dx.z),
    ] {
        close(actual, expected);
    }
    let spin = next
        .principal_frame()
        .properties()
        .angular_velocity_body_checked(next.state().angular_momentum_body())
        .unwrap();
    close(spin.z, 0.2 * (-DT / tau_rot).exp());
    let terminal_slip = terminal.sub(flow());
    let path_integral = tau * -(-DT / tau).exp_m1();
    let translation_heat = c
        * (DT * terminal_slip.norm_squared()
            + 2.0 * path_integral * terminal_slip.dot(deviation)
            + 0.5 * tau * -(-2.0 * DT / tau).exp_m1() * deviation.norm_squared());
    let rotation_heat = 0.5 * inertia * 0.2_f64.powi(2) * -(-2.0 * DT / tau_rot).exp_m1();
    close(report.viscous_heat_j, translation_heat + rotation_heat);
    assert!(report.viscous_heat_j > 0.0);
    assert_eq!(report.thermal.external_heat_j, 0.0);
    assert_eq!(next.thermal().energy_j(), original.thermal().energy_j());
    assert_eq!(next.configuration().inertia().mass(), mass);
    close(report.combined_residual_j, 0.0);
    close(next.thermal().time_s(), DT);
}

#[test]
fn expanding_navier_sphere_reuses_exact_motion_then_preserves_world_momenta() {
    use fs_mbd_material::creeping_flow::navier_slip::NavierSlipBoundary;
    let alpha = 0.02;
    let rho = 0.002;
    let mu = 0.1;
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
    let sphere = CreepingSphere::new(&reservoir, RADIUS, gravity(), flow(), limits())
        .unwrap()
        .with_navier_slip(NavierSlipBoundary::new(&interface).unwrap())
        .unwrap();
    let (next, report) = original
        .advance_creeping_sphere_exact(&sphere, exact_controls(), 4.0 * DT, config(), || Ok(()))
        .unwrap();
    let mass = source.inertia().mass();
    let volume = |r: f64| 4.0 * std::f64::consts::PI * r.powi(3) / 3.0;
    let inertia = |r: f64| 0.4 * mass * r * r;
    let mut velocity = fs_mbd::Vec3::new(0.01, -0.02, 0.03);
    let mut position = fs_mbd::Vec3::new(0.0, 0.0, -0.2);
    let mut spin = 0.2 * inertia(RADIUS);
    let mut parameter_work = 0.0;
    for step in 0..4 {
        let radius = RADIUS * (1.0 + alpha * RATE * f64::from(step) * DT);
        let next_radius = RADIUS * (1.0 + alpha * RATE * f64::from(step + 1) * DT);
        let b_over_r = RADIUS / radius;
        let c = 6.0 * std::f64::consts::PI * mu * radius * (1.0 + 2.0 * b_over_r)
            / (1.0 + 3.0 * b_over_r);
        let c_rot = 8.0 * std::f64::consts::PI * mu * radius.powi(3) / (1.0 + 3.0 * b_over_r);
        let terminal = flow().add(
            gravity()
                .acceleration_world()
                .scale((mass - rho * volume(radius)) / c),
        );
        let deviation = velocity.sub(terminal);
        position = position
            .add(terminal.scale(DT))
            .add(deviation.scale(mass / c * -(-DT * c / mass).exp_m1()));
        velocity = terminal.add(deviation.scale((-DT * c / mass).exp()));
        spin *= (-DT * c_rot / inertia(radius)).exp();
        parameter_work += 0.5 * spin * spin * (1.0 / inertia(next_radius) - 1.0 / inertia(radius))
            + rho
                * (volume(next_radius) - volume(radius))
                * gravity().acceleration_world().dot(position);
    }
    let actual_velocity = next
        .state()
        .center_of_mass_velocity_world(next.principal_frame().properties())
        .unwrap();
    let actual_position = next.state().pose().position_world();
    let actual_spin = next
        .state()
        .pose()
        .orientation()
        .rotate_body_to_world(next.state().angular_momentum_body());
    for (actual, expected) in [
        (actual_velocity.x, velocity.x),
        (actual_velocity.y, velocity.y),
        (actual_velocity.z, velocity.z),
        (actual_position.x, position.x),
        (actual_position.y, position.y),
        (actual_position.z, position.z),
        (actual_spin.z, spin),
    ] {
        close(actual, expected);
    }
    close(report.parameter_work_j, parameter_work);
    close(
        report.radius_after_m,
        RADIUS * (1.0 + alpha * RATE * 4.0 * DT),
    );
    close(report.combined_residual_j, 0.0);
    assert_eq!(next.configuration().inertia().mass(), mass);
    assert_eq!(sphere.navier_slip().unwrap().slip_length_m(), RADIUS);
    assert!(core::ptr::eq(
        sphere.navier_slip().unwrap().source(),
        &interface
    ));
    let (middle, _) = original
        .advance_creeping_sphere_exact(&sphere, exact_controls(), 2.0 * DT, config(), || Ok(()))
        .unwrap();
    let (restarted, _) = middle
        .advance_creeping_sphere_exact(&sphere, exact_controls(), 4.0 * DT, config(), || Ok(()))
        .unwrap();
    assert_eq!(restarted.state(), next.state());
    assert_eq!(restarted.configuration(), next.configuration());
}

#[test]
fn exact_composition_preserves_heat_history_and_late_refusals_publish_nothing() {
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
    let input = PrescribedHeat {
        watts: source.inertia().mass() * CP * RATE,
        time_s: 0.0,
        supplied_j: 0.0,
        refuse_after_s: None,
    };
    let (expected, _) = direct
        .advance_creeping_sphere_exact(&sphere, exact_controls(), 4.0 * DT, config(), || Ok(()))
        .unwrap();
    let (actual, generated, report) = coupled
        .advance_creeping_sphere_exact_with(
            &input,
            &sphere,
            exact_controls(),
            4.0 * DT,
            config(),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(actual.state(), expected.state());
    close(generated.time_s, actual.thermal().time_s());
    close(generated.supplied_j, report.thermal.external_heat_j);
    close(report.combined_residual_j, 0.0);
    let failure = PrescribedHeat {
        refuse_after_s: Some(DT),
        ..input.clone()
    };
    assert!(
        coupled
            .advance_creeping_sphere_exact_with(
                &failure,
                &sphere,
                exact_controls(),
                4.0 * DT,
                config(),
                || Ok(())
            )
            .is_err()
    );
    assert_eq!(failure.time_s, 0.0);
    assert_eq!(failure.supplied_j, 0.0);
    assert_eq!(coupled.thermal().time_s(), 0.0);
    assert_eq!(coupled.state(), direct.state());
    let lift = CreepingSphere::new(&reservoir, RADIUS, gravity(), flow(), limits())
        .unwrap()
        .with_rotational_lift();
    assert!(
        coupled
            .advance_creeping_sphere_exact(&lift, exact_controls(), DT, config(), || Ok(()))
            .is_err()
    );
    assert!(
        coupled
            .advance_creeping_sphere_exact(&sphere, exact_controls(), DT, config(), || Err(
                DepositionError::InvalidInput("checkpoint source refusal")
            ))
            .is_err()
    );
    assert_eq!(coupled.thermal().time_s(), 0.0);
}


fn drive(first: NativeVec3, last: NativeVec3, timescale_s: f64) -> HeatedCreepingAffineForce {
    HeatedCreepingAffineForce {
        force_before_world_n: first,
        force_after_world_n: last,
        minimum_force_variation_timescale_s: timescale_s,
    }
}
fn interpolate_force(first: NativeVec3, last: NativeVec3, fraction: f64) -> NativeVec3 {
    first.scale(1.0 - fraction).add(last.scale(fraction))
}
fn speed(value: NativeVec3) -> f64 {
    value.x.hypot(value.y).hypot(value.z)
}
fn integrate_power(mut f: impl FnMut(f64) -> f64, duration: f64) -> f64 {
    let intervals = 2048;
    let dt = duration / f64::from(intervals);
    let mut sum = f(0.0) + f(duration);
    let mut correction = 0.0;
    for i in 1..intervals {
        let term = f(f64::from(i) * dt) * if i % 2 == 0 { 2.0 } else { 4.0 };
        let adjusted = term - correction;
        let next = sum + adjusted;
        correction = (next - sum) - adjusted;
        sum = next;
    }
    sum * dt / 3.0
}

#[test]
fn heated_affine_force_uses_current_navier_radius_and_true_velocity_power_path() {
    use fs_mbd_material::creeping_flow::navier_slip::NavierSlipBoundary;
    let alpha = 0.02;
    let rho = 0.002;
    let mu = 0.1;
    let duration = 4.0 * DT;
    let first = NativeVec3::new(0.01, -0.005, 0.005);
    let last = NativeVec3::new(-0.005, 0.01, -0.01);
    let material = card(CP, alpha, 500.0);
    let source = reference_body(&material, VolumeMoments::sphere(RADIUS).unwrap());
    let chart = caloric(&material, 1e5, SolidLiquidPhase::Solid);
    let original = body(&source, &material, &chart, RATE, controls(), gravity(), twist());
    let reservoir = fluid(rho, mu);
    let interface = navier_slip_source::constant(material.id(), RADIUS, 400.0);
    let sphere = CreepingSphere::new(&reservoir, RADIUS, gravity(), flow(), limits())
        .unwrap()
        .with_navier_slip(NavierSlipBoundary::new(&interface).unwrap())
        .unwrap();
    let (next, report) = original.advance_creeping_sphere_exact_affine(
        &sphere, drive(first, last, duration), exact_controls(), duration, config(), || Ok(()),
    ).unwrap();

    // Independent affine/exponential ODE solution at each actual thermal epoch,
    // with ordinary scalar quadrature of external power and viscous power.
    let mass = source.inertia().mass();
    let volume = |r: f64| 4.0 * std::f64::consts::PI * r.powi(3) / 3.0;
    let inertia = |r: f64| 0.4 * mass * r * r;
    let mut velocity = NativeVec3::new(0.01, -0.02, 0.03);
    let mut position = NativeVec3::new(0.0, 0.0, -0.2);
    let mut spin = 0.2 * inertia(RADIUS);
    let mut work = 0.0;
    let mut heat = 0.0;
    let mut background = 0.0;
    let mut parameter = 0.0;
    let mut maximum_reynolds: f64 = 0.0;
    for index in 0..4 {
        let time = f64::from(index) * DT;
        let radius = RADIUS * (1.0 + alpha * RATE * time);
        let next_radius = RADIUS * (1.0 + alpha * RATE * (time + DT));
        let slip = RADIUS / radius;
        let c = 6.0 * std::f64::consts::PI * mu * radius
            * (1.0 + 2.0 * slip) / (1.0 + 3.0 * slip);
        let rotary = 8.0 * std::f64::consts::PI * mu * radius.powi(3) / (1.0 + 3.0 * slip);
        let rate = c / mass;
        let before = interpolate_force(first, last, time / duration);
        let after = interpolate_force(first, last, (time + DT) / duration);
        let weight = gravity().acceleration_world().scale(mass - rho * volume(radius));
        let slope = after.sub(before).scale(1.0 / (c * DT));
        let offset = before.add(weight).scale(1.0 / c).sub(slope.scale(1.0 / rate));
        let amplitude = velocity.sub(flow()).sub(offset);
        let path = |t: f64| flow().add(offset).add(slope.scale(t))
            .add(amplitude.scale((-rate * t).exp()));
        work += integrate_power(
            |t| interpolate_force(before, after, t / DT).dot(path(t)), DT,
        );
        heat += integrate_power(|t| c * path(t).sub(flow()).norm_squared(), DT)
            + 0.5 * spin * spin / inertia(radius)
                * -(-2.0 * DT * rotary / inertia(radius)).exp_m1();
        background += integrate_power(|t| -c * path(t).sub(flow()).dot(flow()), DT);
        let acceleration = speed(before.add(weight)).max(speed(after.add(weight))) / mass;
        let initial_speed = speed(velocity.sub(flow()));
        let bound = initial_speed.max(
            (-rate * DT).exp() * initial_speed
                + acceleration * -(-rate * DT).exp_m1() / rate,
        );
        maximum_reynolds = maximum_reynolds.max(2.0 * radius * rho * bound / mu);
        position = position.add(flow().add(offset).scale(DT))
            .add(slope.scale(0.5 * DT * DT))
            .add(amplitude.scale(-(-rate * DT).exp_m1() / rate));
        velocity = path(DT);
        spin *= (-DT * rotary / inertia(radius)).exp();
        parameter += 0.5 * spin * spin * (1.0 / inertia(next_radius) - 1.0 / inertia(radius))
            + rho * (volume(next_radius) - volume(radius))
                * gravity().acceleration_world().dot(position);
        maximum_reynolds = maximum_reynolds
            .max(2.0 * next_radius * rho * speed(velocity.sub(flow())) / mu);
    }
    let actual_velocity = next.state()
        .center_of_mass_velocity_world(next.principal_frame().properties()).unwrap();
    let actual_position = next.state().pose().position_world();
    let actual_spin = next.state().pose().orientation()
        .rotate_body_to_world(next.state().angular_momentum_body());
    for (actual, expected) in [
        (actual_velocity.x, velocity.x), (actual_velocity.y, velocity.y),
        (actual_velocity.z, velocity.z), (actual_position.x, position.x),
        (actual_position.y, position.y), (actual_position.z, position.z),
        (actual_spin.z, spin),
    ] { close(actual, expected); }
    assert!(work.abs() > 1e-9 && heat > 0.0 && background.abs() > 1e-9);
    close(report.applied_work_j / work, 1.0);
    close(report.viscous_heat_j / heat, 1.0);
    close(report.background_work_j / background, 1.0);
    close(report.parameter_work_j, parameter);
    close(report.maximum_translation_reynolds_bound, maximum_reynolds);
    assert!(report.maximum_translation_reynolds_bound >
        report.load_before.translation_reynolds.max(report.load_after.translation_reynolds));
    let final_radius = RADIUS * (1.0 + alpha * RATE * duration);
    close(report.maximum_drive_diffusion_ratio, final_radius.powi(2) * rho / (mu * duration));
    close(report.radius_after_m, final_radius);
    close(report.thermal.external_heat_j, mass * CP * RATE * duration);
    close(report.combined_residual_j, 0.0);
    assert_eq!(next.configuration().inertia().mass(), mass);
    assert_eq!(sphere.navier_slip().unwrap().slip_length_m(), RADIUS);
    assert!(core::ptr::eq(sphere.navier_slip().unwrap().source(), &interface));
}

#[test]
fn affine_history_survives_thermal_partition_and_true_endpoint_continuation() {
    let material = card(CP, 0.0, 500.0);
    let source = reference_body(&material, VolumeMoments::sphere(RADIUS).unwrap());
    let chart = caloric(&material, 1e5, SolidLiquidPhase::Solid);
    let mut coarse_controls = controls();
    coarse_controls.maximum_step_s = 2.0 * DT;
    let original = body(&source, &material, &chart, RATE, controls(), gravity(), twist());
    let coarse = body(&source, &material, &chart, RATE, coarse_controls, gravity(), twist());
    let reservoir = fluid(0.002, 0.1);
    let sphere = CreepingSphere::new(&reservoir, RADIUS, gravity(), flow(), limits()).unwrap();
    let duration = 4.0 * DT;
    let first = NativeVec3::new(0.002, -0.001, 0.001);
    let last = NativeVec3::new(-0.001, 0.002, -0.002);
    let (whole, audit) = original.advance_creeping_sphere_exact_affine(
        &sphere, drive(first, last, duration), exact_controls(), duration, config(), || Ok(()),
    ).unwrap();
    let (partitioned, partition_audit) = coarse.advance_creeping_sphere_exact_affine(
        &sphere, drive(first, last, duration), exact_controls(), duration, config(), || Ok(()),
    ).unwrap();
    let middle_force = interpolate_force(first, last, 0.5);
    let (middle, first_audit) = original.advance_creeping_sphere_exact_affine(
        &sphere, drive(first, middle_force, duration), exact_controls(),
        2.0 * DT, config(), || Ok(()),
    ).unwrap();
    let (restarted, last_audit) = middle.advance_creeping_sphere_exact_affine(
        &sphere, drive(middle_force, last, duration), exact_controls(), duration, config(), || Ok(()),
    ).unwrap();
    for candidate in [&partitioned, &restarted] {
        let a = candidate.state();
        let b = whole.state();
        for (actual, expected) in [
            (a.linear_momentum_world().x, b.linear_momentum_world().x),
            (a.linear_momentum_world().y, b.linear_momentum_world().y),
            (a.linear_momentum_world().z, b.linear_momentum_world().z),
            (a.pose().position_world().x, b.pose().position_world().x),
            (a.pose().position_world().y, b.pose().position_world().y),
            (a.pose().position_world().z, b.pose().position_world().z),
        ] { close(actual, expected); }
    }
    close(partition_audit.applied_work_j / audit.applied_work_j, 1.0);
    close((first_audit.applied_work_j + last_audit.applied_work_j) / audit.applied_work_j, 1.0);
    close((first_audit.viscous_heat_j + last_audit.viscous_heat_j) / audit.viscous_heat_j, 1.0);
    assert_eq!(partition_audit.maximum_drive_diffusion_ratio, audit.maximum_drive_diffusion_ratio);
    assert_eq!(last_audit.maximum_drive_diffusion_ratio, audit.maximum_drive_diffusion_ratio);

    let (unforced, old_audit) = original.advance_creeping_sphere_exact(
        &sphere, exact_controls(), duration, config(), || Ok(()),
    ).unwrap();
    let (zero, zero_audit) = original.advance_creeping_sphere_exact_affine(
        &sphere, drive(NativeVec3::ZERO, NativeVec3::ZERO, 0.0),
        exact_controls(), duration, config(), || Ok(()),
    ).unwrap();
    assert_eq!(zero.state(), unforced.state());
    assert_eq!(zero_audit, old_audit);
    assert_eq!(zero_audit.applied_work_j, 0.0);
}

#[test]
fn driven_heat_observer_composes_once_and_rolls_back_after_late_failure() {
    let material = card(CP, 0.02, 500.0);
    let source = reference_body(&material, VolumeMoments::sphere(RADIUS).unwrap());
    let chart = caloric(&material, 1e5, SolidLiquidPhase::Solid);
    let direct = body(&source, &material, &chart, RATE, controls(), gravity(), twist());
    let coupled = body(&source, &material, &chart, 0.0, controls(), gravity(), twist());
    let reservoir = fluid(0.002, 0.1);
    let sphere = CreepingSphere::new(&reservoir, RADIUS, gravity(), flow(), limits()).unwrap();
    let force = drive(
        NativeVec3::new(0.002, -0.001, 0.001),
        NativeVec3::new(-0.001, 0.002, -0.002),
        4.0 * DT,
    );
    let input = PrescribedHeat {
        watts: source.inertia().mass() * CP * RATE,
        time_s: 0.0,
        supplied_j: 0.0,
        refuse_after_s: None,
    };
    let (expected, expected_audit) = direct.advance_creeping_sphere_exact_affine(
        &sphere, force, exact_controls(), 4.0 * DT, config(), || Ok(()),
    ).unwrap();
    let (actual, generated, audit) = coupled.advance_creeping_sphere_exact_affine_with(
        &input, &sphere, force, exact_controls(), 4.0 * DT, config(), || Ok(()),
    ).unwrap();
    assert_eq!(actual.state(), expected.state());
    assert_eq!(actual.configuration(), expected.configuration());
    close(generated.supplied_j, audit.thermal.external_heat_j);
    close(audit.applied_work_j, expected_audit.applied_work_j);
    close(audit.combined_residual_j, 0.0);
    let failure = PrescribedHeat { refuse_after_s: Some(DT), ..input.clone() };
    assert!(coupled.advance_creeping_sphere_exact_affine_with(
        &failure, &sphere, force, exact_controls(), 4.0 * DT, config(), || Ok(()),
    ).is_err());
    assert_eq!(failure.time_s, 0.0);
    assert_eq!(failure.supplied_j, 0.0);
    assert_eq!(coupled.thermal().time_s(), 0.0);
    assert_eq!(coupled.state(), direct.state());
}

#[test]
fn expansion_rechecks_the_physical_drive_timescale_at_updated_endpoints() {
    let material = card(CP, 0.02, 500.0);
    let source = reference_body(&material, VolumeMoments::sphere(RADIUS).unwrap());
    let chart = caloric(&material, 1e5, SolidLiquidPhase::Solid);
    let original = body(&source, &material, &chart, RATE, controls(), gravity(), twist());
    let rho = 0.002;
    let mu = 0.001;
    let reservoir = fluid(rho, mu);
    let sphere = CreepingSphere::new(&reservoir, RADIUS, gravity(), flow(), limits()).unwrap();
    let physical_timescale = RADIUS.powi(2) * rho / (mu * 0.048);
    let first = NativeVec3::new(1e-4, 0.0, 0.0);
    let last = NativeVec3::new(1.1e-4, 0.0, 0.0);
    let (_, prefix) = original.advance_creeping_sphere_exact_affine(
        &sphere, drive(first, interpolate_force(first, last, 0.5), physical_timescale),
        exact_controls(), 2.0 * DT, config(), || Ok(()),
    ).unwrap();
    close(prefix.maximum_drive_diffusion_ratio, 0.048 * 1.02_f64.powi(2));
    assert!(prefix.maximum_diffusion_ratio < limits().maximum_diffusion_ratio);
    assert!(matches!(original.advance_creeping_sphere_exact_affine(
        &sphere, drive(first, last, physical_timescale),
        exact_controls(), 4.0 * DT, config(), || Ok(()),
    ), Err(HeatedMotionError::CreepingFlow(CreepingFlowError::UnsteadyDrive { .. }))));
    assert_eq!(original.thermal().time_s(), 0.0);
    assert_eq!(original.configuration().inertia().mass(), source.inertia().mass());
}
