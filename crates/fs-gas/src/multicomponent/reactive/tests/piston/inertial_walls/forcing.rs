//! G1/G3/G4: chemical storage, body work and two finite-wall reactions.

use super::*;
use crate::body_force::{BodyForce, ForceControl};

const FORCE_CONTROL: ForceControl = ForceControl { maximum_mach: 80.0 };

#[test]
fn g1_reacting_accelerated_gas_closes_independent_complete_wall_balances() {
    let mechanism = mechanism(8.0, 30.0, 0.0);
    let mut flow = closed(&mechanism, 8, chemistry());
    let before = flow.flow().totals();
    let initial_walls = [0.2, 0.5].map(|mass| Piston::new(mass, 0.0).unwrap());
    let mut walls = initial_walls;
    let mut loads = wall_loads(1.0);
    loads[0].external_force_n = 0.003;
    loads[1].external_force_n = -0.002;
    let acceleration = 3.0;
    let step = flow
        .step_with_inertial_walls_and_body_force(
            &mut walls,
            Time::new(0.006),
            ReactiveStepControl::default(),
            loads,
            wall_controls(),
            BodyForce::UniformAcceleration(acceleration),
            FORCE_CONTROL,
            || false,
        )
        .unwrap();
    let after = flow.flow().totals();
    let h = step.flow.elapsed.value();
    let reference_change: f64 = (0..2)
        .map(|k| {
            (after.species_masses[k] - before.species_masses[k])
                * flow.flow().gas().energy_references()[k]
        })
        .sum();
    let physical_change = after.energy - before.energy + reference_change;
    let kinetic: f64 = walls
        .iter()
        .map(|wall| 0.5 * wall.mass_kg() * wall.velocity_m_s().powi(2))
        .sum();
    let positions = flow.flow().mesh().positions();
    let displacement = [positions[0], positions[8] - 1.0];
    let exterior_work = AREA * (displacement[1] - displacement[0]);
    let load_work: f64 = (0..2)
        .map(|i| loads[i].external_force_n * displacement[i])
        .sum();
    let body = step.flow.transport.source_change;
    close(body.impulse, before.mass() * acceleration * h, 3e-11);
    assert!(body.work > 0.0);
    assert_eq!(body, step.walls.gas.source_change);
    close(step.gas_energy_change_j, physical_change, 3e-9);
    close(
        physical_change + kinetic + exterior_work - load_work - body.work,
        0.0,
        3e-8,
    );
    close(
        step.total_energy_residual_j,
        0.0,
        wall_controls().energy_tolerance_j,
    );
    let combined_momentum = after.momentum - before.momentum
        + walls
            .iter()
            .map(|wall| wall.mass_kg() * wall.velocity_m_s())
            .sum::<f64>();
    let external_impulse = h * (loads[0].external_force_n + loads[1].external_force_n);
    close(
        combined_momentum - external_impulse - body.impulse,
        0.0,
        3e-11,
    );
    close(
        step.combined_momentum_residual_n_s,
        0.0,
        wall_controls().momentum_tolerance_n_s,
    );
    close(after.mass(), before.mass(), 3e-12);
    assert!(step.flow.chemistry.species_change_kg[1] > 0.0);
    assert!(step.flow.chemistry.referenced_energy_change_j > 0.0);
    assert_eq!(step.flow.elapsed, step.walls.gas.elapsed);
    assert_eq!(flow.flow().time(), step.flow.elapsed);
}

#[test]
fn g1_isomerizing_free_fall_keeps_gas_and_walls_on_the_same_accelerating_motion() {
    let mechanism = mechanism(0.0, 2.0, 0.0);
    let mut flow = closed(&mechanism, 8, chemistry());
    let acceleration = -2.0;
    let mut walls = [0.2, 0.5].map(|mass| Piston::new(mass, 0.0).unwrap());
    let loads = walls.map(|wall| InertialWallLoad {
        exterior_pressure_pa: 1.0,
        external_force_n: wall.mass_kg() * acceleration,
    });
    let step = flow
        .step_with_inertial_walls_and_body_force(
            &mut walls,
            Time::new(0.001),
            ReactiveStepControl::default(),
            loads,
            wall_controls(),
            BodyForce::UniformAcceleration(acceleration),
            FORCE_CONTROL,
            || false,
        )
        .unwrap();
    let h = step.flow.elapsed.value();
    for wall in walls {
        close(wall.velocity_m_s(), acceleration * h, 3e-10);
    }
    close(
        flow.flow().mesh().positions()[0],
        0.5 * acceleration * h * h,
        3e-12,
    );
    close(
        flow.flow().mesh().positions()[8],
        1.0 + 0.5 * acceleration * h * h,
        3e-12,
    );
    let expected_a = 0.5 + 0.5 * (-4.0 * h).exp();
    for i in 0..8 {
        let q = flow.flow().cell(i).unwrap();
        let p = flow.flow().gas().primitive(q, || false).unwrap();
        close(p.velocity, acceleration * h, 3e-9);
        close(p.temperature, 1.0, 3e-9);
        close(p.density, 1.0, 3e-10);
        close(q.partial_densities[0], expected_a, 2e-6);
    }
    close(
        step.total_energy_residual_j,
        0.0,
        wall_controls().energy_tolerance_j,
    );
}

#[test]
fn g4_forced_reacting_wall_refusal_preserves_geometry_chemistry_and_both_velocities() {
    let mechanism = mechanism(8.0, 30.0, 0.0);
    let h = Time::new(0.001);
    let loads = wall_loads(1.0);
    let initial_walls = [0.2, 0.5].map(|mass| Piston::new(mass, 0.0).unwrap());
    let mut ordinary = closed(&mechanism, 8, chemistry());
    let mut zero = closed(&mechanism, 8, chemistry());
    let mut a = initial_walls;
    let mut b = initial_walls;
    let expected = ordinary
        .step_with_inertial_walls(
            &mut a,
            h,
            ReactiveStepControl::default(),
            loads,
            wall_controls(),
            || false,
        )
        .unwrap();
    let actual = zero
        .step_with_inertial_walls_and_body_force(
            &mut b,
            h,
            ReactiveStepControl::default(),
            loads,
            wall_controls(),
            BodyForce::UniformAcceleration(0.0),
            FORCE_CONTROL,
            || false,
        )
        .unwrap();
    assert_eq!(actual, expected);
    assert_eq!(zero.accepted.states, ordinary.accepted.states);
    assert_eq!(a, b);

    let force = BodyForce::UniformAcceleration(3.0);
    let mut reference = closed(&mechanism, 8, chemistry());
    let mut reference_walls = initial_walls;
    let mut polls = 0;
    let expected = reference
        .step_with_inertial_walls_and_body_force(
            &mut reference_walls,
            h,
            ReactiveStepControl::default(),
            loads,
            wall_controls(),
            force,
            FORCE_CONTROL,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut flow = closed(&mechanism, 8, chemistry());
    let mut walls = initial_walls;
    let before = flow.accepted.states.clone();
    let positions = flow.flow().mesh().positions().to_vec();
    let volumes = flow.flow().mesh().volumes().to_vec();
    let mut attempted = 0;
    assert_eq!(
        flow.step_with_inertial_walls_and_body_force(
            &mut walls,
            h,
            ReactiveStepControl::default(),
            loads,
            wall_controls(),
            force,
            FORCE_CONTROL,
            || {
                attempted += 1;
                attempted >= polls - 1
            },
        )
        .unwrap_err(),
        ReactiveError::Cancelled
    );
    assert_eq!(flow.accepted.states, before);
    assert_eq!(flow.flow().mesh().positions(), positions);
    assert_eq!(flow.flow().mesh().volumes(), volumes);
    assert_eq!(flow.flow().time().value(), 0.0);
    assert_eq!(walls, initial_walls);
    let replay = flow
        .step_with_inertial_walls_and_body_force(
            &mut walls,
            h,
            ReactiveStepControl::default(),
            loads,
            wall_controls(),
            force,
            FORCE_CONTROL,
            || false,
        )
        .unwrap();
    assert_eq!(replay.flow.elapsed, expected.flow.elapsed);
    for (actual, expected) in flow.accepted.states.iter().zip(&reference.accepted.states) {
        close(*actual, *expected, 3e-9);
    }
    for (actual, expected) in walls.iter().zip(reference_walls) {
        close(actual.velocity_m_s(), expected.velocity_m_s(), 3e-9);
    }
}
