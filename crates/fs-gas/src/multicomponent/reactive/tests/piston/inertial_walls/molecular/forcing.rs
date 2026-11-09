//! G1/G3/G4: actual body work with viscous, reacting, thermally active walls.

use super::*;
use crate::body_force::{BodyForce, ForceControl};

const FORCE_CONTROL: ForceControl = ForceControl { maximum_mach: 80.0 };

fn forced_heat(
    flow: &mut ReactingEuler1d,
    operator: &mut DissipationOperator<ConstantTransport>,
    thermal: &mut PistonThermalEnds<'_>,
    walls: &mut [Piston; 2],
    force: BodyForce<'_>,
    cancelled: impl FnMut() -> bool,
) -> Result<
    crate::dynamics::inertial_walls::InertialWallsStep<ThermalReactiveStepReport>,
    ReactiveError,
> {
    flow.step_with_inertial_walls_dissipation_heat_and_body_force(
        operator,
        thermal,
        walls,
        Time::new(0.005),
        ReactiveStepControl::default(),
        DissipationControl::default(),
        WallControl::default(),
        wall_loads(1.0),
        wall_controls(),
        force,
        FORCE_CONTROL,
        cancelled,
    )
}

#[test]
fn g1_acceleration_closes_actual_chemical_thermal_and_both_wall_storages() {
    for acceleration in [-3.0, 2.0] {
        let (mut flow, mut operator, mut thermal, mut walls) = thermal_fixture();
        let before = flow.flow().totals();
        let contacts = *thermal.contacts();
        let initial_walls = walls;
        let initial_positions = flow.flow().mesh().positions().to_vec();
        let report = forced_heat(
            &mut flow,
            &mut operator,
            &mut thermal,
            &mut walls,
            BodyForce::UniformAcceleration(acceleration),
            || false,
        )
        .unwrap();
        let after = flow.flow().totals();
        let h = report.gas.reactive.elapsed.value();
        let body = report.gas.reactive.transport.source_change;
        close(body.impulse, before.mass() * acceleration * h, 3e-11);
        let store_change = 0.4 * (temperature(thermal.contacts()[0]) - temperature(contacts[0]))
            + 0.7 * (temperature(thermal.contacts()[1]) - temperature(contacts[1]));
        let kinetic: f64 = (0..2)
            .map(|i| {
                0.5 * walls[i].mass_kg()
                    * (walls[i].velocity_m_s().powi(2) - initial_walls[i].velocity_m_s().powi(2))
            })
            .sum();
        let right = flow.flow().cell_count();
        let displacement =
            [0, right].map(|i| flow.flow().mesh().positions()[i] - initial_positions[i]);
        let exterior_work = AREA * (displacement[1] - displacement[0]);
        let gas_change = physical_change(flow.flow().gas(), before, after);
        close(report.gas.gas_energy_change_j, gas_change, 4e-10);
        close(
            gas_change + store_change + kinetic + exterior_work - body.work,
            0.0,
            4e-8,
        );
        close(
            after.momentum - before.momentum
                + (0..2)
                    .map(|i| {
                        walls[i].mass_kg()
                            * (walls[i].velocity_m_s() - initial_walls[i].velocity_m_s())
                    })
                    .sum::<f64>()
                - body.impulse,
            0.0,
            4e-11,
        );
        // The wall reaction contains native pressure and longitudinal stress,
        // never the distributed external acceleration or conductive heat.
        let molecular = report.gas.reactive.dissipation.unwrap();
        let pressure = report.gas.reactive.transport;
        let changes = [
            [
                pressure.left_boundary_change,
                pressure.right_boundary_change,
            ],
            [
                molecular.first.left_boundary_change,
                molecular.first.right_boundary_change,
            ],
            [
                molecular.second.left_boundary_change,
                molecular.second.right_boundary_change,
            ],
        ];
        for end in 0..2 {
            close(
                report.gas_impulse_n_s[end],
                -changes.iter().map(|pair| pair[end].momentum).sum::<f64>(),
                3e-14,
            );
            for pair in &changes {
                close(
                    pair[end].energy,
                    report.midpoint_velocity_m_s[end] * pair[end].momentum,
                    3e-13,
                );
            }
        }
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
        assert!(report.gas.reactive.chemistry.species_change_kg[1] > 0.0);
        assert!(report.gas.heat.heat_to_gas_j.abs() > 1e-8);
        assert_eq!(report.gas.heat.external_heat_to_gas_j, 0.0);
        assert_eq!(operator.mesh(), flow.flow().mesh());
        assert_eq!(thermal.mesh(), flow.flow().mesh());
    }
}

#[test]
fn g1_forced_molecular_trial_replays_native_stages_and_actual_moving_volumes() {
    let mechanism = mechanism(40.0, 500.0, 0.0);
    let mut flow = closed(&mechanism, 4, chemistry());
    let mut replay = closed(&mechanism, 4, chemistry());
    let mut operator = molecular(flow.flow());
    let mut replay_operator = molecular(replay.flow());
    let mut walls = [2.0, 3.0].map(|mass| Piston::new(mass, 0.0).unwrap());
    let force_density = [1.0, 2.0, -0.5, 3.0];
    let report = flow
        .step_with_inertial_walls_dissipation_and_body_force(
            &mut operator,
            &mut walls,
            Time::new(0.2),
            ReactiveStepControl::default(),
            DissipationControl::default(),
            wall_loads(1.0),
            wall_controls(),
            BodyForce::ForceDensity(&force_density),
            FORCE_CONTROL,
            || false,
        )
        .unwrap();
    assert!(report.gas.elapsed.value() < 0.2);
    let h = report.gas.elapsed.value();
    let speeds = report.midpoint_velocity_m_s;
    replay_operator
        .admit_motion(&replay.accepted, &mut || false)
        .unwrap();
    replay_operator
        .advance_moving_trial(
            &mut replay.accepted,
            Time::new(0.5 * h),
            DissipationControl::default(),
            speeds,
            &mut || false,
        )
        .unwrap();
    let mut chemistry_report = ChemistryReport::default();
    react_cells(
        &mut replay.accepted,
        &mut replay.reactor,
        &mut replay.amounts,
        0.5 * h,
        &mut chemistry_report,
        &mut || false,
    )
    .unwrap();
    let old_volumes = replay.flow().mesh().volumes().to_vec();
    let transport = replay
        .accepted
        .step_with_moving_walls_and_body_force(
            Time::new(h),
            ReactiveStepControl::default().transport,
            speeds,
            BodyForce::ForceDensity(&force_density),
            FORCE_CONTROL,
            || false,
        )
        .unwrap();
    assert_eq!(transport.elapsed.value(), h);
    let impulse: f64 = force_density
        .iter()
        .enumerate()
        .map(|(i, f)| 0.5 * h * f * (old_volumes[i] + replay.flow().mesh().volumes()[i]))
        .sum();
    close(transport.source_change.impulse, impulse, 3e-12);
    assert_ne!(replay.flow().mesh().volumes(), old_volumes.as_slice());
    react_cells(
        &mut replay.accepted,
        &mut replay.reactor,
        &mut replay.amounts,
        0.5 * h,
        &mut chemistry_report,
        &mut || false,
    )
    .unwrap();
    replay_operator
        .advance_moving_trial(
            &mut replay.accepted,
            Time::new(0.5 * h),
            DissipationControl::default(),
            speeds,
            &mut || false,
        )
        .unwrap();
    assert_eq!(flow.flow().mesh(), replay.flow().mesh());
    assert_eq!(flow.flow().time(), replay.flow().time());
    assert_eq!(report.gas.transport.source_change, transport.source_change);
    for (&actual, &expected) in flow.accepted.states.iter().zip(&replay.accepted.states) {
        close(actual, expected, 3e-12);
    }
}

#[test]
fn g3_zero_body_force_preserves_native_molecular_heat_path_bitwise() {
    let (mut ordinary, mut da, mut ta, mut wa) = thermal_fixture();
    let (mut zero, mut db, mut tb, mut wb) = thermal_fixture();
    for _ in 0..2 {
        let expected = heat_step(&mut ordinary, &mut da, &mut ta, &mut wa, || false).unwrap();
        let actual = forced_heat(
            &mut zero,
            &mut db,
            &mut tb,
            &mut wb,
            BodyForce::Acceleration(&[0.0; 4]),
            || false,
        )
        .unwrap();
        assert_eq!(actual, expected);
        assert_eq!(zero.accepted.states, ordinary.accepted.states);
        assert_eq!(zero.flow().mesh(), ordinary.flow().mesh());
        assert_eq!(wb, wa);
        assert_eq!(tb.contacts(), ta.contacts());
    }
}

#[test]
fn g4_forced_heat_late_cancellation_and_bad_source_preserve_all_owners() {
    let force = BodyForce::UniformAcceleration(3.0);
    let (mut reference, mut dr, mut tr, mut wr) = thermal_fixture();
    forced_heat(&mut reference, &mut dr, &mut tr, &mut wr, force, || false).unwrap();
    let states = reference.accepted.states.clone();
    let mesh = reference.flow().mesh().clone();
    let time = reference.flow().time();
    let initial_walls = wr;
    let contacts = *tr.contacts();
    let last_heat = tr.last_report();
    let mut polls = 0;
    let expected = forced_heat(&mut reference, &mut dr, &mut tr, &mut wr, force, || {
        polls += 1;
        false
    })
    .unwrap();
    for cutoff in [polls / 2, polls] {
        let (mut flow, mut operator, mut thermal, mut walls) = thermal_fixture();
        forced_heat(
            &mut flow,
            &mut operator,
            &mut thermal,
            &mut walls,
            force,
            || false,
        )
        .unwrap();
        let mut seen = 0;
        let error = forced_heat(
            &mut flow,
            &mut operator,
            &mut thermal,
            &mut walls,
            force,
            || {
                seen += 1;
                seen == cutoff
            },
        )
        .unwrap_err();
        assert_eq!(error, ReactiveError::Cancelled);
        assert_eq!(seen, cutoff);
        assert_eq!(flow.accepted.states, states);
        assert_eq!(flow.flow().time(), time);
        assert_eq!(flow.flow().mesh(), &mesh);
        assert_eq!(operator.mesh(), &mesh);
        assert_eq!(thermal.mesh(), &mesh);
        assert_eq!(walls, initial_walls);
        assert_eq!(*thermal.contacts(), contacts);
        assert_eq!(thermal.last_report(), last_heat);
        assert!(
            forced_heat(
                &mut flow,
                &mut operator,
                &mut thermal,
                &mut walls,
                BodyForce::ForceDensity(&[1.0]),
                || false,
            )
            .is_err()
        );
        assert_eq!(flow.accepted.states, states);
        assert_eq!(walls, initial_walls);
        assert_eq!(*thermal.contacts(), contacts);
        let resumed = forced_heat(
            &mut flow,
            &mut operator,
            &mut thermal,
            &mut walls,
            force,
            || false,
        )
        .unwrap();
        assert_eq!(resumed, expected);
        assert_eq!(flow.accepted.states, reference.accepted.states);
        assert_eq!(walls, wr);
        assert_eq!(thermal.contacts(), tr.contacts());
    }
}
