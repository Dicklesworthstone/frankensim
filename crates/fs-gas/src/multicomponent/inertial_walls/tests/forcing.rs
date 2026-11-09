//! G1/G3/G4: native pressure/wall reactions and separate distributed acceleration.

use super::*;

const FORCE_CONTROL: ForceControl = ForceControl { maximum_mach: 80.0 };

fn uniform(velocity: f64) -> MixtureEuler1d {
    let gas = field([0.0; 3]).gas().clone();
    let states = (0..24)
        .map(|_| gas.conserved(&[0.3, 0.7, 0.0], velocity, 2.0).unwrap())
        .collect();
    MixtureEuler1d::new(gas, Length::new(1.0 / 24.0), AREA, states, CLOSED, || false).unwrap()
}

#[test]
fn g1_free_fall_of_gas_and_both_walls_retains_exact_uniform_translation() {
    for acceleration in [-3.0, 3.0] {
        let speed = 0.2;
        let mut flow = uniform(speed);
        let mass = flow.totals().mass();
        let mut walls = [0.25, 1.1].map(|mass| Piston::new(mass, speed).unwrap());
        let before_momentum = momentum(&flow, walls);
        let before_energy = energy(&flow, walls);
        let pressure = flow
            .gas()
            .primitive(flow.cell(0).unwrap(), || false)
            .unwrap()
            .pressure;
        let loads = walls.map(|wall| InertialWallLoad {
            exterior_pressure_pa: pressure,
            external_force_n: wall.mass_kg() * acceleration,
        });
        let report = flow
            .step_with_inertial_walls_and_body_force(
                &mut walls,
                Time::new(0.0025),
                StepControl::default(),
                loads,
                controls(),
                BodyForce::UniformAcceleration(acceleration),
                FORCE_CONTROL,
                || false,
            )
            .unwrap();
        let h = report.gas.elapsed.value();
        let displacement = h * (speed + 0.5 * acceleration * h);
        let final_speed = speed + acceleration * h;
        close(
            report.gas.source_change.impulse,
            mass * acceleration * h,
            3e-12,
        );
        close(
            report.gas.source_change.work,
            mass * acceleration * displacement,
            3e-12,
        );
        for wall in walls {
            close(wall.velocity_m_s(), final_speed, 3e-11);
        }
        close(flow.mesh().positions()[0], displacement, 3e-13);
        close(flow.mesh().positions()[24], 1.0 + displacement, 3e-13);
        for i in 0..24 {
            let cell = flow.cell(i).unwrap();
            let p = flow.gas().primitive(cell, || false).unwrap();
            close(p.density, 1.0, 3e-11);
            close(p.temperature, 2.0, 3e-10);
            close(p.velocity, final_speed, 3e-10);
            assert_eq!(cell.partial_densities[2], 0.0);
        }
        let total_mass = mass + walls.iter().map(|w| w.mass_kg()).sum::<f64>();
        close(
            momentum(&flow, walls) - before_momentum,
            total_mass * acceleration * h,
            3e-11,
        );
        close(
            energy(&flow, walls) - before_energy,
            total_mass * acceleration * displacement,
            3e-10,
        );
        close(
            report.combined_momentum_residual_n_s,
            0.0,
            controls().momentum_tolerance_n_s,
        );
        close(
            report.total_energy_residual_j,
            0.0,
            controls().energy_tolerance_j,
        );
    }
}

#[test]
fn g1_force_density_uses_actual_moving_volumes_and_native_pressure_proposals() {
    let mut flow = field([0.0; 3]);
    let mut replay = field([0.0; 3]);
    let mut walls = [0.25, 1.1].map(|mass| Piston::new(mass, 0.0).unwrap());
    let before_volume: f64 = flow.mesh().volumes().iter().sum();
    let force = BodyForce::UniformForceDensity(2.0);
    let report = flow
        .step_with_inertial_walls_and_body_force(
            &mut walls,
            Time::new(0.0025),
            StepControl::default(),
            [InertialWallLoad::default(); 2],
            controls(),
            force,
            FORCE_CONTROL,
            || false,
        )
        .unwrap();
    let after_volume: f64 = flow.mesh().volumes().iter().sum();
    close(
        report.gas.source_change.impulse,
        report.gas.elapsed.value() * (before_volume + after_volume),
        3e-12,
    );
    assert!(after_volume > before_volume);
    let reference = replay
        .step_with_moving_walls_and_body_force(
            report.gas.elapsed,
            StepControl::default(),
            report.midpoint_velocity_m_s,
            force,
            FORCE_CONTROL,
            || false,
        )
        .unwrap();
    for (actual, expected) in flow.states.iter().zip(&replay.states) {
        close(*actual, *expected, 3e-11);
    }
    close(
        report.gas.source_change.work,
        reference.source_change.work,
        3e-12,
    );
    for i in 0..2 {
        let gas_side = [
            reference.left_boundary_change,
            reference.right_boundary_change,
        ][i];
        close(report.gas_impulse_n_s[i], -gas_side.momentum, 3e-12);
        close(
            report.kinetic_energy_change_j[i] + gas_side.energy,
            0.0,
            3e-11,
        );
    }
    close(
        report.combined_momentum_residual_n_s,
        0.0,
        controls().momentum_tolerance_n_s,
    );
    close(
        report.total_energy_residual_j,
        0.0,
        controls().energy_tolerance_j,
    );
}

#[test]
fn g4_zero_forcing_and_late_cancel_preserve_all_native_wall_state() {
    let h = Time::new(0.0025);
    let loads = [InertialWallLoad::default(); 2];
    let initial_walls = [0.25, 1.1].map(|mass| Piston::new(mass, 0.0).unwrap());
    let mut ordinary = field([0.0; 3]);
    let mut zero = field([0.0; 3]);
    let mut a = initial_walls;
    let mut b = initial_walls;
    let expected = ordinary
        .step_with_inertial_walls(&mut a, h, StepControl::default(), loads, controls(), || {
            false
        })
        .unwrap();
    let actual = zero
        .step_with_inertial_walls_and_body_force(
            &mut b,
            h,
            StepControl::default(),
            loads,
            controls(),
            BodyForce::Acceleration(&[0.0; 24]),
            FORCE_CONTROL,
            || false,
        )
        .unwrap();
    assert_eq!(actual, expected);
    assert_eq!(a, b);
    assert_eq!(zero.states, ordinary.states);
    assert!(zero.body_force_stage.is_none());

    let mut reference = field([0.0; 3]);
    let mut reference_walls = initial_walls;
    let force = BodyForce::UniformAcceleration(3.0);
    let mut polls = 0;
    let accepted = reference
        .step_with_inertial_walls_and_body_force(
            &mut reference_walls,
            h,
            StepControl::default(),
            loads,
            controls(),
            force,
            FORCE_CONTROL,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut flow = field([0.0; 3]);
    let mut walls = initial_walls;
    let original = flow.states.clone();
    let positions = flow.mesh().positions().to_vec();
    let volumes = flow.mesh().volumes().to_vec();
    let mut attempted = 0;
    let result = flow.step_with_inertial_walls_and_body_force(
        &mut walls,
        h,
        StepControl::default(),
        loads,
        controls(),
        force,
        FORCE_CONTROL,
        || {
            attempted += 1;
            attempted >= polls - 1
        },
    );
    assert!(matches!(result, Err(PistonError::Cancelled)));
    assert_eq!(flow.states, original);
    assert_eq!(flow.mesh().positions(), positions);
    assert_eq!(flow.mesh().volumes(), volumes);
    assert_eq!(flow.time().value(), 0.0);
    assert_eq!(walls, initial_walls);
    assert!(
        flow.step_with_inertial_walls_and_body_force(
            &mut walls,
            h,
            StepControl::default(),
            loads,
            controls(),
            BodyForce::Acceleration(&[1.0]),
            FORCE_CONTROL,
            || false,
        )
        .is_err()
    );
    assert_eq!(flow.states, original);
    let retry = flow
        .step_with_inertial_walls_and_body_force(
            &mut walls,
            h,
            StepControl::default(),
            loads,
            controls(),
            force,
            FORCE_CONTROL,
            || false,
        )
        .unwrap();
    assert_eq!(retry.gas.elapsed, accepted.gas.elapsed);
    for (actual, expected) in flow.states.iter().zip(&reference.states) {
        close(*actual, *expected, 3e-10);
    }
}
