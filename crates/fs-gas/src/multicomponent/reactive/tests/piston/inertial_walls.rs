//! G1/G3/G4: original reaction energy driving both finite end walls.
use super::*;
use crate::dynamics::inertial_walls::{InertialWallLoad, InertialWallsControl};

mod molecular;

fn wall_controls() -> InertialWallsControl {
    let control = piston_controls();
    InertialWallsControl {
        midpoint_velocity_bracket_m_s: [control.midpoint_velocity_bracket_m_s; 2],
        momentum_tolerance_n_s: control.momentum_tolerance_n_s,
        energy_tolerance_j: control.energy_tolerance_j,
        maximum_iterations: control.maximum_iterations,
    }
}

fn wall_loads(pressure: f64) -> [InertialWallLoad; 2] {
    [InertialWallLoad {
        exterior_pressure_pa: pressure,
        external_force_n: 0.0,
    }; 2]
}

#[test]
fn g1_exothermic_two_wall_recoil_closes_physical_energy_momentum_and_source_gauge() {
    let mut reference = None;
    for datum in [0.0, 128.0] {
        let mechanism = mechanism(8.0, 30.0, datum);
        let mut flow = closed(&mechanism, 8, chemistry());
        let before = flow.flow().totals();
        let mut walls = [0.2, 0.5].map(|mass| Piston::new(mass, 0.0).unwrap());
        let report = flow
            .step_with_inertial_walls(
                &mut walls,
                Time::new(0.01),
                ReactiveStepControl::default(),
                wall_loads(1.0),
                wall_controls(),
                || false,
            )
            .unwrap();
        let after = flow.flow().totals();
        assert_eq!(report.flow.elapsed, report.walls.gas.elapsed);
        assert_eq!(flow.flow().time(), report.flow.elapsed);
        assert_eq!(report.flow.chemistry.cell_solves, 16);
        assert!(walls[0].velocity_m_s() < 0.0);
        assert!(walls[1].velocity_m_s() > 0.0);
        assert!(
            after.momentum.abs() > 1e-12,
            "unequal walls leave resolved gas momentum"
        );
        assert!(after.species_masses[0] < before.species_masses[0]);
        assert!(after.species_masses[1] > 0.0);
        let positions = flow.flow().mesh().positions();
        assert!(positions[0] < 0.0);
        assert!(positions[8] > 1.0);
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
        // Both exterior pressures are 1 Pa; their signed impulses cancel.
        // Their work opposes expansion at BOTH represented mesh endpoints.
        let exterior_work = AREA * ((positions[8] - 1.0) - positions[0]);
        close(physical_change + kinetic + exterior_work, 0.0, 2e-8);
        close(report.gas_energy_change_j, physical_change, 2e-8);
        close(
            after.momentum - before.momentum
                + walls
                    .iter()
                    .map(|wall| wall.mass_kg() * wall.velocity_m_s())
                    .sum::<f64>(),
            0.0,
            3e-11,
        );
        close(
            report.total_energy_residual_j,
            0.0,
            wall_controls().energy_tolerance_j,
        );
        close(
            report.combined_momentum_residual_n_s,
            0.0,
            wall_controls().momentum_tolerance_n_s,
        );
        close(after.mass(), before.mass(), 3e-12);
        close(
            flow.flow().mesh().volumes().iter().sum(),
            AREA * (positions[8] - positions[0]),
            5e-14,
        );
        for k in 0..2 {
            assert_eq!(report.flow.transport.boundary_change.species_masses[k], 0.0);
            close(
                after.species_masses[k] - before.species_masses[k],
                report.flow.chemistry.species_change_kg[k],
                3e-12,
            );
            close(
                report.flow.conservation_residual.species_masses[k],
                0.0,
                3e-12,
            );
        }
        let physical: Vec<_> = (0..8)
            .map(|i| {
                let p = flow
                    .flow()
                    .gas()
                    .primitive(flow.flow().cell(i).unwrap(), || false)
                    .unwrap();
                [p.temperature, p.velocity, p.pressure]
            })
            .collect();
        let value = (walls.map(|wall| wall.velocity_m_s()), physical);
        if let Some((velocities, states)) = &reference {
            for (a, b) in value.0.into_iter().zip(*velocities) {
                close(a, b, 2e-9);
            }
            for (a, b) in value.1.iter().zip(states) {
                for (&a, &b) in a.iter().zip(b) {
                    close(a, b, 2e-9);
                }
            }
        } else {
            reference = Some(value);
        }
    }
}

#[test]
fn g1_two_wall_acoustic_restart_replays_chemistry_on_actual_moved_volumes() {
    let mechanism = mechanism(40.0, 500.0, 0.0);
    let mut flow = closed(&mechanism, 4, chemistry());
    let initial_walls = [2.0, 3.0].map(|mass| Piston::new(mass, 0.0).unwrap());
    let mut walls = initial_walls;
    let report = flow
        .step_with_inertial_walls(
            &mut walls,
            Time::new(0.2),
            ReactiveStepControl::default(),
            wall_loads(1.0),
            wall_controls(),
            || false,
        )
        .unwrap();
    assert!(
        report.flow.coupling_retries > 0,
        "reaction heat must shorten the acoustic duration"
    );
    let dt = report.flow.elapsed.value();
    assert_eq!(report.flow.elapsed, report.walls.gas.elapsed);
    assert_eq!(report.flow.chemistry.cell_solves, 8);
    let mut replay = closed(&mechanism, 4, chemistry());
    let original_volumes = replay.flow().mesh().volumes().to_vec();
    let mut replay_walls = initial_walls;
    let mut sources = ChemistryReport::default();
    react_cells(
        &mut replay.accepted,
        &mut replay.reactor,
        &mut replay.amounts,
        0.5 * dt,
        &mut sources,
        &mut || false,
    )
    .unwrap();
    let middle = replay
        .accepted
        .step_with_inertial_walls(
            &mut replay_walls,
            Time::new(dt),
            StepControl::default(),
            wall_loads(1.0),
            wall_controls(),
            || false,
        )
        .unwrap();
    assert_eq!(middle.gas.elapsed, report.flow.elapsed);
    assert!(middle.displacement_m[0] < 0.0);
    assert!(middle.displacement_m[1] > 0.0);
    assert_ne!(replay.flow().mesh().volumes(), original_volumes.as_slice());
    react_cells(
        &mut replay.accepted,
        &mut replay.reactor,
        &mut replay.amounts,
        0.5 * dt,
        &mut sources,
        &mut || false,
    )
    .unwrap();
    for i in 0..2 {
        close(
            walls[i].velocity_m_s(),
            replay_walls[i].velocity_m_s(),
            2e-10,
        );
    }
    for (&a, &b) in flow.accepted.states.iter().zip(&replay.accepted.states) {
        close(a, b, 2e-10);
    }
    for (&a, &b) in flow
        .flow()
        .mesh()
        .positions()
        .iter()
        .zip(replay.flow().mesh().positions())
    {
        close(a, b, 2e-12);
    }
    for k in 0..2 {
        close(
            report.flow.chemistry.species_change_kg[k],
            sources.species_change_kg[k],
            2e-10,
        );
    }
    close(
        report.total_energy_residual_j,
        0.0,
        wall_controls().energy_tolerance_j,
    );
}

#[test]
fn g4_second_chemistry_refusal_rolls_back_both_compressing_walls_and_mesh() {
    let mechanism = mechanism(0.0, 2.0, 0.0);
    let chemistry = IsochoricKineticsSolve {
        maximum_temperature_k: 1.0000001,
        ..chemistry()
    };
    let duration = Time::new(0.001);
    let original_walls = [0.1, 0.15].map(|mass| Piston::new(mass, 0.0).unwrap());
    let mut stages = closed(&mechanism, 8, chemistry);
    let mut stage_walls = original_walls;
    let mut sources = ChemistryReport::default();
    react_cells(
        &mut stages.accepted,
        &mut stages.reactor,
        &mut stages.amounts,
        0.5 * duration.value(),
        &mut sources,
        &mut || false,
    )
    .unwrap();
    let middle = stages
        .accepted
        .step_with_inertial_walls(
            &mut stage_walls,
            duration,
            StepControl::default(),
            wall_loads(2.0),
            wall_controls(),
            || false,
        )
        .unwrap();
    assert_eq!(middle.gas.elapsed, duration);
    assert!(middle.displacement_m[0] > 0.0);
    assert!(middle.displacement_m[1] < 0.0);
    assert_eq!(
        react_cells(
            &mut stages.accepted,
            &mut stages.reactor,
            &mut stages.amounts,
            0.5 * duration.value(),
            &mut sources,
            &mut || false,
        ),
        Err(ReactiveError::Invalid(
            "cell temperature is outside chemistry interval"
        )),
    );
    let mut flow = closed(&mechanism, 8, chemistry);
    let before = flow.accepted.states.clone();
    let mesh = flow.flow().mesh().clone();
    let mut walls = original_walls;
    let result = flow.step_with_inertial_walls(
        &mut walls,
        duration,
        ReactiveStepControl {
            maximum_coupling_retries: 0,
            ..ReactiveStepControl::default()
        },
        wall_loads(2.0),
        wall_controls(),
        || false,
    );
    assert!(matches!(result, Err(ReactiveError::RetryLimit { last, .. })
        if *last == ReactiveError::Invalid("cell temperature is outside chemistry interval")));
    assert_eq!(flow.accepted.states, before);
    assert_eq!(flow.flow().mesh(), &mesh);
    assert_eq!(flow.flow().time(), Time::new(0.0));
    assert_eq!(walls, original_walls);
}

#[test]
fn g4_cancellation_after_two_wall_motion_preserves_all_owners_and_retry() {
    let mechanism = mechanism(10.0, 30.0, 0.0);
    let duration = Time::new(0.002);
    let initial_walls = [0.2, 0.5].map(|mass| Piston::new(mass, 0.0).unwrap());
    let mut reference = closed(&mechanism, 8, chemistry());
    let mut reference_walls = initial_walls;
    reference
        .step_with_inertial_walls(
            &mut reference_walls,
            duration,
            ReactiveStepControl::default(),
            wall_loads(1.0),
            wall_controls(),
            || false,
        )
        .unwrap();
    let before = reference.accepted.states.clone();
    let mesh = reference.flow().mesh().clone();
    let time = reference.flow().time();
    let original_walls = reference_walls;
    let mut polls = 0;
    let expected = reference
        .step_with_inertial_walls(
            &mut reference_walls,
            duration,
            ReactiveStepControl::default(),
            wall_loads(1.0),
            wall_controls(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    for cutoff in [1, polls / 2, polls] {
        let mut flow = closed(&mechanism, 8, chemistry());
        let mut walls = initial_walls;
        flow.step_with_inertial_walls(
            &mut walls,
            duration,
            ReactiveStepControl::default(),
            wall_loads(1.0),
            wall_controls(),
            || false,
        )
        .unwrap();
        let mut seen = 0;
        assert_eq!(
            flow.step_with_inertial_walls(
                &mut walls,
                duration,
                ReactiveStepControl::default(),
                wall_loads(1.0),
                wall_controls(),
                || {
                    seen += 1;
                    seen == cutoff
                },
            ),
            Err(ReactiveError::Cancelled),
        );
        assert_eq!(seen, cutoff);
        assert_eq!(flow.accepted.states, before);
        assert_eq!(flow.flow().mesh(), &mesh);
        assert_eq!(flow.flow().time(), time);
        assert_eq!(walls, original_walls);
        let recovered = flow
            .step_with_inertial_walls(
                &mut walls,
                duration,
                ReactiveStepControl::default(),
                wall_loads(1.0),
                wall_controls(),
                || false,
            )
            .unwrap();
        assert_eq!(recovered.flow.elapsed, expected.flow.elapsed);
        for i in 0..2 {
            close(
                walls[i].velocity_m_s(),
                reference_walls[i].velocity_m_s(),
                2e-10,
            );
        }
        for (&a, &b) in flow.accepted.states.iter().zip(&reference.accepted.states) {
            close(a, b, 2e-10);
        }
        for (&a, &b) in flow
            .flow()
            .mesh()
            .positions()
            .iter()
            .zip(reference.flow().mesh().positions())
        {
            close(a, b, 2e-13);
        }
    }
}
