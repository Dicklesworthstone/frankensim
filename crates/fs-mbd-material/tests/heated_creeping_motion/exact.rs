//! Analytic frozen-radius motion on the existing accepted thermal timeline.
use super::*;
use fs_mbd_material::thermal::heating::motion::creeping_flow::HeatedExactCreepingControls;

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
