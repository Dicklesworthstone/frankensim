//! G1/G4: external drive, source feedback and finite heat share one acceptance.
use super::*;
use fs_couple::friction_heat::directional::driven::*;
use fs_mbd::Wrench;
use fs_tribo::tangential::elliptic::driven::{EllipticDrivenControls, EllipticDrivenMode};
use fs_tribo::tangential::elliptic::relaxation::EllipticRelaxationControls;

fn driven_controls() -> EllipticDrivenControls {
    EllipticDrivenControls {
        maximum_step_s: 0.002,
        minimum_step_s: 1e-10,
        maximum_substeps: 100_000,
        rest_velocity_tolerance_m_s: 1e-10,
        solve: controls().mechanical.contact,
        relaxation: EllipticRelaxationControls {
            time_tolerance_s: 1e-11,
            travel_tolerance_m: 1e-11,
            maximum_evaluations: 40_000,
            maximum_root_iterations: 96,
        },
    }
}

fn resting_mechanics(
    common_speed: f64,
) -> (
    [RigidBodyState; 3],
    [MassProperties; 3],
    RigidDirectionalContact,
) {
    let (mut bodies, properties, contact) = mechanics(false);
    for i in 0..2 {
        bodies[i] = RigidBodyState::new(
            bodies[i].pose(),
            Vec3::new((2.0 + i as f64) * common_speed, 0.0, 0.0),
            Vec3::ZERO,
        )
        .unwrap();
    }
    (bodies, properties, contact)
}

fn drive(force: f64) -> [Wrench; 2] {
    [
        Wrench {
            force_world: Vec3::new(force, 0.0, 0.0),
            torque_body: Vec3::ZERO,
        },
        Wrench::ZERO,
    ]
}

#[test]
fn static_hold_supplies_common_body_work_without_generating_friction_heat() {
    let (a, ca) = material(0);
    let (b, cb) = material(1);
    let cards = [a, b];
    let curves = [ca, cb];
    let q = query();
    let card = interface(&cards, [scalar(0.2), scalar(0.4)], 400.0);
    let law = law(&card, &q, [1.0, 0.0]);
    let mut history = law.initial_state();
    let mut mode = EllipticDrivenMode::Resting;
    let mut heat = thermal(&curves, [300.0; 2]);
    let initial_cells = heat.cells().clone();
    let route = FrictionHeatRoute::from_link(&heat, 7, CONDUCTIVITY).unwrap();
    let (mut bodies, props, contact) = resting_mechanics(1.0);
    let initial = bodies;
    let dt = 0.01;
    let traction = 120.0 / 101.0; // a/W = (2N/2kg)/(101/120 kg^-1).
    let work = 2.0 * dt + 0.5 * (2.0 - traction) * dt * dt;
    let report = apply_source_directional_driven_heat(
        &mut bodies,
        &props,
        &mut history,
        &mut mode,
        &law,
        contact,
        drive(2.0),
        &mut heat,
        &route,
        binding(&cards, &curves, &card, &q),
        dt,
        controls(),
        driven_controls(),
        &[0.0; 2],
        300.0,
        |_, p| transport(&cards, p),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(mode, EllipticDrivenMode::Resting);
    assert_eq!(report.mechanical.driven.velocity_after_m_s, [0.0; 2]);
    close(report.mechanical.slip_velocity_after_m_s[0], 0.0);
    close(
        report.mechanical.driven.traction_impulse_ns[0],
        traction * dt,
    );
    close(
        bodies[0].linear_momentum_world().x,
        2.0 + (2.0 - traction) * dt,
    );
    close(bodies[1].linear_momentum_world().x, 3.0 + traction * dt);
    for body in &bodies[..2] {
        close(body.angular_momentum_body().y, -0.1 * traction * dt);
    }
    close(report.external_mechanical_work_j, work);
    close(report.mechanical.kinetic_energy_change_j, work);
    close(report.mechanical.driven.drive_work_j, 0.0);
    close(report.heat.dissipated_j, 0.0);
    close(report.actual_distance_m, 0.0);
    assert_eq!(report.net_slip_m, [0.0; 2]);
    assert_eq!(heat.cells(), &initial_cells);
    close(heat.time_s(), dt);
    close(report.combined_energy_residual_j, 0.0);
    for source in history.friction_sources() {
        close(source.state_coordinate("speed").unwrap(), 0.0);
        close(source.state_coordinate("path").unwrap(), 0.25);
    }
    assert_eq!(bodies[2], initial[2]);
    assert_eq!(bodies[0].pose(), initial[0].pose());
    assert_eq!(bodies[1].pose(), initial[1].pose());
}

#[test]
fn breakaway_heat_and_external_work_match_independent_constant_force_mechanics() {
    let (a, ca) = material(0);
    let (b, cb) = material(1);
    let cards = [a, b];
    let curves = [ca, cb];
    let q = query();
    let card = interface(&cards, [scalar(0.2), scalar(0.4)], 400.0);
    let law = law(&card, &q, [1.0, 0.0]);
    let mut history = law.initial_state();
    let mut mode = EllipticDrivenMode::Resting;
    let mut heat = thermal(&curves, [300.0; 2]);
    let route = FrictionHeatRoute::from_link(&heat, 7, CONDUCTIVITY).unwrap();
    let (mut bodies, props, contact) = resting_mechanics(0.0);
    let initial = bodies;
    let dt = 0.01;
    let mobility = 101.0 / 120.0;
    let acceleration = 10.0 - mobility * 2.0;
    let distance = 0.5 * acceleration * dt * dt;
    let friction_heat = 2.0 * distance;
    // External force acts at A's centre: x_A = (20-2)dt^2/(2*m_A).
    let external_work = 20.0 * 0.5 * 9.0 * dt * dt;
    let report = apply_source_directional_driven_heat(
        &mut bodies,
        &props,
        &mut history,
        &mut mode,
        &law,
        contact,
        drive(20.0),
        &mut heat,
        &route,
        binding(&cards, &curves, &card, &q),
        dt,
        controls(),
        driven_controls(),
        &[0.0; 2],
        300.0,
        |_, p| transport(&cards, p),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(mode, EllipticDrivenMode::Sliding);
    close(report.mechanical.mobility_per_kg[0][0], mobility);
    close(report.mechanical.driving_acceleration_m_s2[0], 10.0);
    close(
        report.mechanical.slip_velocity_after_m_s[0],
        acceleration * dt,
    );
    close(report.mechanical.driven.traction_impulse_ns[0], 2.0 * dt);
    close(bodies[0].linear_momentum_world().x, 18.0 * dt);
    close(bodies[1].linear_momentum_world().x, 2.0 * dt);
    for body in &bodies[..2] {
        close(body.angular_momentum_body().y, -0.2 * dt);
    }
    close(report.actual_distance_m, distance);
    close(report.net_slip_m[0], distance);
    close(report.net_slip_m[1], 0.0);
    close(report.heat.dissipated_j, friction_heat);
    close(report.external_mechanical_work_j, external_work);
    close(
        report.mechanical.kinetic_energy_change_j,
        external_work - friction_heat,
    );
    let actual_kinetic: f64 = (0..2)
        .map(|i| bodies[i].kinetic_energy(props[i]).unwrap())
        .sum();
    let actual_thermal: f64 = (0..2)
        .map(|i| CAPACITY[i] * (heat.temperature_k(CELLS[i]).unwrap() - 300.0))
        .sum();
    close(actual_kinetic + actual_thermal, external_work);
    close(actual_thermal, friction_heat);
    close(report.heat.body_heat_j[0], 0.25 * friction_heat);
    close(report.heat.body_heat_j[1], 0.75 * friction_heat);
    close(
        report.heat.interface_temperature_k,
        300.0 + friction_heat / 0.04,
    );
    close(report.pressure_pa, 10.0);
    close(report.combined_energy_residual_j, 0.0);
    for source in history.friction_sources() {
        close(source.state_coordinate("speed").unwrap(), distance / dt);
        close(source.state_coordinate("path").unwrap(), 0.25 + distance);
        close(source.state_coordinate("pressure").unwrap(), 10.0);
    }
    assert_eq!(bodies[2], initial[2]);
}

#[test]
fn driven_temperature_speed_and_path_feedback_retain_original_claims() {
    let (a, ca) = material(0);
    let (b, cb) = material(1);
    let cards = [a, b];
    let curves = [ca, cb];
    let q = query();
    let dt = 0.01;
    let mobility = 101.0 / 120.0;
    for (case, kinetic) in [
        curve("Ti", Temperature::DIMS, [(299.0, 0.3), (301.0, 0.1)]),
        curve("speed", SPEED, [(0.0, 0.1), (20.0, 0.7)]),
        curve("path", LENGTH, [(0.0, 0.2), (10.0, 0.6)]),
    ]
    .into_iter()
    .enumerate()
    {
        let card = interface(&cards, [kinetic, scalar(0.4)], 400.0);
        let law = law(&card, &q, [1.0, 0.0]);
        let mut history = law.initial_state();
        let initial = history.clone();
        let mut mode = EllipticDrivenMode::Resting;
        let mut heat = thermal(&curves, [300.0; 2]);
        let route = FrictionHeatRoute::from_link(&heat, 7, CONDUCTIVITY).unwrap();
        let (mut bodies, props, contact) = resting_mechanics(0.0);
        // Eliminate both finite cells and their massless face independently.
        // r = dt^2(a-WF)/2, D=F*r, Ti=300+D/.04, speed=r/dt.
        let mut lower = 0.0;
        let mut upper = 4.0;
        for _ in 0..80 {
            let force = 0.5 * (lower + upper);
            let path = 0.5 * dt * dt * (10.0 - mobility * force);
            let expected_force = match case {
                0 => 2.0 - force * path / 0.04,
                1 => 1.0 + 0.3 * path / dt,
                _ => 2.0 + 0.4 * (0.25 + path),
            };
            if force > expected_force {
                upper = force;
            } else {
                lower = force;
            }
        }
        let force = 0.5 * (lower + upper);
        let path = 0.5 * dt * dt * (10.0 - mobility * force);
        let dissipated = force * path;
        let report = apply_source_directional_driven_heat(
            &mut bodies,
            &props,
            &mut history,
            &mut mode,
            &law,
            contact,
            drive(20.0),
            &mut heat,
            &route,
            binding(&cards, &curves, &card, &q),
            dt,
            controls(),
            driven_controls(),
            &[0.0; 2],
            300.0,
            |_, p| transport(&cards, p),
            &CancelGate::new(),
        )
        .unwrap();
        assert_eq!(mode, EllipticDrivenMode::Sliding);
        close(report.mechanical.driven.traction_impulse_ns[0] / dt, force);
        close(report.actual_distance_m, path);
        close(report.heat.dissipated_j, dissipated);
        close(
            report.heat.interface_temperature_k,
            300.0 + dissipated / 0.04,
        );
        close(
            report.external_mechanical_work_j,
            5.0 * (20.0 - force) * dt * dt,
        );
        close(report.combined_energy_residual_j, 0.0);
        assert!(report.source_motion_residual_m <= 0.0);
        assert!(report.source_motion_residual_m.abs() < 2e-12);
        close(
            history.friction_sources()[0]
                .property("kinetic-x")
                .unwrap()
                .value_si(),
            force / 10.0,
        );
        for (i, source) in history.friction_sources().iter().enumerate() {
            close(source.state_coordinate("speed").unwrap(), path / dt);
            close(source.state_coordinate("path").unwrap(), 0.25 + path);
            close(
                source.state_coordinate("Ti").unwrap(),
                300.0 + dissipated / 0.04,
            );
            close(
                source.state_coordinate("Ta").unwrap(),
                300.0 + dissipated / 0.08,
            );
            close(
                source.state_coordinate("Tb").unwrap(),
                300.0 + dissipated / 0.08,
            );
            close(
                report.source[i].query().axes()["T"],
                300.0 + dissipated / 0.08,
            );
            for property in source.properties() {
                assert_eq!(
                    property.answer().receipt.selected,
                    initial.friction_sources()[i]
                        .property(property.requirement().name())
                        .unwrap()
                        .answer()
                        .receipt
                        .selected
                );
                card.claims()
                    .verify_receipt(&property.answer().receipt)
                    .unwrap();
            }
        }
    }
}

#[test]
fn driven_failures_preserve_mode_bodies_and_heat_then_retry_and_zero_load_work() {
    let (a, ca) = material(0);
    let (b, cb) = material(1);
    let cards = [a, b];
    let curves = [ca, cb];
    let q = query();
    let card = interface(&cards, [scalar(0.2), scalar(0.4)], 400.0);
    let law = law(&card, &q, [1.0, 0.0]);
    let initial = law.initial_state();
    let base = thermal(&curves, [300.0; 2]);
    let route = FrictionHeatRoute::from_link(&base, 7, CONDUCTIVITY).unwrap();
    let (original, props, contact) = resting_mechanics(0.0);
    for failure in 0..5 {
        let mut bodies = original;
        let mut history = initial.clone();
        let mut mode = EllipticDrivenMode::Resting;
        let mut heat = base.clone();
        let gate = CancelGate::new();
        let mut solve = controls();
        let mut selected = contact;
        if failure == 2 {
            solve.maximum_trials = 1;
        }
        if failure == 3 {
            solve.motion.distance_bracket_m = [0.0, 1e-5];
        }
        if failure == 4 {
            selected.nominal_pressure_pa = Some(10.0);
        }
        let result = apply_source_directional_driven_heat(
            &mut bodies,
            &props,
            &mut history,
            &mut mode,
            &law,
            selected,
            drive(20.0),
            &mut heat,
            &route,
            binding(&cards, &curves, &card, &q),
            0.01,
            solve,
            driven_controls(),
            &[0.0; 2],
            300.0,
            |_, p| {
                if failure == 0 {
                    return Err(ThermalMaterialError::Invalid("driven transport refused"));
                }
                if failure == 1 {
                    gate.request();
                }
                transport(&cards, p)
            },
            &gate,
        );
        assert!(
            result.is_err(),
            "driven failure {failure} unexpectedly accepted"
        );
        assert_eq!(bodies, original);
        assert_eq!(history, initial);
        assert_eq!(mode, EllipticDrivenMode::Resting);
        assert_eq!(heat.cells(), base.cells());
        assert_eq!(heat.links(), base.links());
        assert_eq!(heat.time_s(), base.time_s());
    }
    let mut bodies = original;
    let mut history = initial.clone();
    let mut mode = EllipticDrivenMode::Resting;
    let mut heat = base.clone();
    let report = apply_source_directional_driven_heat(
        &mut bodies,
        &props,
        &mut history,
        &mut mode,
        &law,
        contact,
        drive(20.0),
        &mut heat,
        &route,
        binding(&cards, &curves, &card, &q),
        0.01,
        controls(),
        driven_controls(),
        &[0.0; 2],
        300.0,
        |_, p| transport(&cards, p),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(mode, EllipticDrivenMode::Sliding);
    assert_ne!(history, initial);
    assert_ne!(bodies, original);
    close(report.combined_energy_residual_j, 0.0);

    // A zero load means zero sourced pressure and zero friction, but explicit
    // body work and existing hot-to-cold conduction still advance together.
    let mut contact = contact;
    contact.normal_load_n = 0.0;
    let mut bodies = original;
    let mut history = initial;
    let mut mode = EllipticDrivenMode::Resting;
    let mut heat = thermal(&curves, [280.0, 320.0]);
    let report = apply_source_directional_driven_heat(
        &mut bodies,
        &props,
        &mut history,
        &mut mode,
        &law,
        contact,
        drive(20.0),
        &mut heat,
        &route,
        binding(&cards, &curves, &card, &q),
        0.01,
        controls(),
        driven_controls(),
        &[0.0; 2],
        310.0,
        |_, p| transport(&cards, p),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(mode, EllipticDrivenMode::Sliding);
    close(report.pressure_pa, 0.0);
    close(report.heat.dissipated_j, 0.0);
    close(report.external_mechanical_work_j, 0.01);
    close(report.mechanical.kinetic_energy_change_j, 0.01);
    close(report.actual_distance_m, 0.0005);
    close(report.heat.body_heat_j[0], 0.3);
    close(report.heat.body_heat_j[1], -0.3);
    close(heat.temperature_k(CELLS[0]).unwrap(), 295.0);
    close(heat.temperature_k(CELLS[1]).unwrap(), 315.0);
    close(report.combined_energy_residual_j, 0.0);
    for source in history.friction_sources() {
        close(source.state_coordinate("pressure").unwrap(), 0.0);
    }
    assert_eq!(bodies[2], original[2]);
}

#[path = "driven/march.rs"]
mod march;
