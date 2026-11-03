//! Independent mechanical/thermal balances and moving-owner rollback.
use super::super::molecular::molecular;
use super::*;
use crate::multicomponent::dissipation::{
    ConstantTransport, DissipationControl, DissipationOperator,
};
use crate::multicomponent::wall_exchange::{WallContact, WallControl, WallDrag, WallThermal};

fn nonuniform(reference: f64) -> ReactingEuler1d {
    let mechanism = mechanism(8.0, 30.0, reference);
    let mut flow = field(&mechanism, 8, 0.0, |x| {
        let a = 0.7 + 0.1 * (2.0 * core::f64::consts::PI * x).sin();
        [a, 1.0 - a]
    });
    flow.boundaries = Boundaries::Ends {
        left: Boundary::Reflecting,
        right: Boundary::Reflecting,
    };
    ReactingEuler1d::new(flow, mechanism, chemistry(), || false).unwrap()
}

fn physical_change(gas: &MixtureGas, before: Totals, after: Totals) -> f64 {
    sum(
        core::iter::once(after.energy - before.energy).chain((0..gas.species_count()).map(|k| {
            (after.species_masses[k] - before.species_masses[k]) * gas.energy_references()[k]
        })),
    )
}

fn finite(temperature_k: f64, capacity: f64, conductance: f64) -> WallContact {
    WallContact {
        drag: WallDrag::None,
        thermal: WallThermal::Lumped {
            temperature_k,
            heat_capacity_j_per_k: capacity,
            conductance_w_per_k: conductance,
        },
    }
}

fn temperature(contact: WallContact) -> f64 {
    match contact.thermal {
        WallThermal::Lumped { temperature_k, .. } => temperature_k,
        _ => panic!("expected a finite thermal store"),
    }
}

fn heat_step(
    flow: &mut ReactingEuler1d,
    operator: &mut DissipationOperator<ConstantTransport>,
    thermal: &mut PistonThermalEnds<'_>,
    walls: &mut [Piston; 2],
    cancelled: impl FnMut() -> bool,
) -> Result<
    crate::dynamics::inertial_walls::InertialWallsStep<ThermalReactiveStepReport>,
    ReactiveError,
> {
    flow.step_with_inertial_walls_dissipation_and_heat(
        operator,
        thermal,
        walls,
        Time::new(0.005),
        ReactiveStepControl::default(),
        DissipationControl::default(),
        WallControl::default(),
        wall_loads(1.0),
        wall_controls(),
        cancelled,
    )
}

#[test]
fn g1_two_viscous_walls_close_actual_energy_momentum_and_source_gauge() {
    let mut a = nonuniform(0.0);
    let mut b = nonuniform(1000.0);
    let mut da = molecular(a.flow());
    let mut db = molecular(b.flow());
    let mut wa = [0.2, 0.5].map(|mass| Piston::new(mass, 0.0).unwrap());
    let mut wb = wa;
    let before = a.flow().totals();
    let mut exterior_work = 0.0;
    let mut viscous_impulse = [0.0; 2];
    for _ in 0..3 {
        let report = a
            .step_with_inertial_walls_and_dissipation(
                &mut da,
                &mut wa,
                Time::new(0.005),
                ReactiveStepControl::default(),
                DissipationControl::default(),
                wall_loads(1.0),
                wall_controls(),
                || false,
            )
            .unwrap();
        let shifted = b
            .step_with_inertial_walls_and_dissipation(
                &mut db,
                &mut wb,
                Time::new(0.005),
                ReactiveStepControl::default(),
                DissipationControl::default(),
                wall_loads(1.0),
                wall_controls(),
                || false,
            )
            .unwrap();
        let halves = report.gas.dissipation.unwrap();
        for half in [halves.first, halves.second] {
            assert_eq!(half.elapsed.value(), 0.5 * report.gas.elapsed.value());
            assert_eq!(half.boundary_change.species_masses[..2], [0.0; 2]);
            for (end, change) in [half.left_boundary_change, half.right_boundary_change]
                .into_iter()
                .enumerate()
            {
                close(
                    change.energy,
                    report.midpoint_velocity_m_s[end] * change.momentum,
                    3e-14,
                );
                viscous_impulse[end] += change.momentum.abs();
            }
        }
        let pressure = [
            report.gas.transport.left_boundary_change.momentum,
            report.gas.transport.right_boundary_change.momentum,
        ];
        let first = [
            halves.first.left_boundary_change.momentum,
            halves.first.right_boundary_change.momentum,
        ];
        let second = [
            halves.second.left_boundary_change.momentum,
            halves.second.right_boundary_change.momentum,
        ];
        for end in 0..2 {
            close(
                report.gas_impulse_n_s[end],
                -pressure[end] - first[end] - second[end],
                3e-14,
            );
            close(
                report.mechanical_energy_residual_j[end],
                0.0,
                wall_controls().energy_tolerance_j,
            );
            close(
                report.momentum_residual_n_s[end],
                0.0,
                wall_controls().momentum_tolerance_n_s,
            );
            close(wa[end].velocity_m_s(), wb[end].velocity_m_s(), 3e-9);
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
        assert_eq!(da.mesh(), a.flow().mesh());
        assert_eq!(db.mesh(), b.flow().mesh());
        close(
            report.gas.elapsed.value(),
            shifted.gas.elapsed.value(),
            1e-15,
        );
        exterior_work += report.exterior_pressure_work_j.into_iter().sum::<f64>();
    }
    let after = a.flow().totals();
    let kinetic: f64 = wa
        .iter()
        .map(|w| 0.5 * w.mass_kg() * w.velocity_m_s().powi(2))
        .sum();
    close(
        physical_change(a.flow().gas(), before, after) + kinetic + exterior_work,
        0.0,
        8e-8,
    );
    // Equal exterior pressures have opposite signed impulses; no fixed support.
    close(
        after.momentum - before.momentum
            + wa.iter()
                .map(|w| w.mass_kg() * w.velocity_m_s())
                .sum::<f64>(),
        0.0,
        2e-10,
    );
    assert!(viscous_impulse.into_iter().all(|impulse| impulse > 1e-10));
    assert!(wa[0].velocity_m_s() < 0.0 && wa[1].velocity_m_s() > 0.0);
    assert!(after.species_masses[0] < before.species_masses[0]);
    for i in 0..a.flow().cell_count() {
        let x = a
            .flow()
            .gas()
            .primitive(a.flow().cell(i).unwrap(), || false)
            .unwrap();
        let y = b
            .flow()
            .gas()
            .primitive(b.flow().cell(i).unwrap(), || false)
            .unwrap();
        for (left, right) in [x.temperature, x.pressure, x.velocity].into_iter().zip([
            y.temperature,
            y.pressure,
            y.velocity,
        ]) {
            close(left, right, 3e-9);
        }
    }
}

#[test]
fn g1_both_finite_end_stores_drive_recoil_with_independent_heat_and_energy_oracles() {
    let mechanism = mechanism(0.0, 0.0, 0.0);
    for wall_t in [3.0, 0.4] {
        let mut flow = closed(&mechanism, 1, chemistry());
        let mut operator = molecular(flow.flow());
        let mut thermal = PistonThermalEnds::new(
            flow.flow(),
            [finite(wall_t, 0.4, 1.0), finite(wall_t, 0.7, 1.5)],
            [None, None],
            || false,
        )
        .unwrap();
        let mut walls = [0.2, 0.5].map(|mass| Piston::new(mass, 0.0).unwrap());
        let before = flow.flow().totals();
        let mut exterior_work = 0.0;
        for step in 0..3 {
            let report =
                heat_step(&mut flow, &mut operator, &mut thermal, &mut walls, || false).unwrap();
            if step == 0 {
                // Independent sequential backward Euler for both finite stores.
                // One gas cell has mass .1 kg and cv=2.5 J/(kg K).
                let half = 0.5 * report.gas.reactive.elapsed.value();
                let q_left = half * (wall_t - 1.0) / (1.0 + half * (1.0 / 0.25 + 1.0 / 0.4));
                let after_left = 1.0 + q_left / 0.25;
                let q_right = half * 1.5 * (wall_t - after_left)
                    / (1.0 + half * 1.5 * (1.0 / 0.25 + 1.0 / 0.7));
                close(report.gas.heat.first.ends[0].heat_to_gas_j, q_left, 2e-10);
                close(report.gas.heat.first.ends[1].heat_to_gas_j, q_right, 2e-10);
            }
            for half in [report.gas.heat.first, report.gas.heat.second] {
                for end in half.ends {
                    assert_eq!(end.wall_impulse_ns, 0.0);
                    assert_eq!(end.drag_heating_j, 0.0);
                }
            }
            assert_eq!(report.gas.heat.external_heat_to_gas_j, 0.0);
            close(report.gas.heat.energy_residual_j, 0.0, 4e-10);
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
            for residual in report.mechanical_energy_residual_j {
                close(residual, 0.0, wall_controls().energy_tolerance_j);
            }
            assert_eq!(thermal.last_report(), Some(report.gas.heat));
            assert_eq!(thermal.mesh(), flow.flow().mesh());
            assert_eq!(operator.mesh(), flow.flow().mesh());
            exterior_work += report.exterior_pressure_work_j.into_iter().sum::<f64>();
        }
        let sign = (wall_t - 1.0).signum();
        assert!(sign * walls[0].velocity_m_s() < 0.0);
        assert!(sign * walls[1].velocity_m_s() > 0.0);
        let store_change = 0.4 * (temperature(thermal.contacts()[0]) - wall_t)
            + 0.7 * (temperature(thermal.contacts()[1]) - wall_t);
        let kinetic: f64 = walls
            .iter()
            .map(|w| 0.5 * w.mass_kg() * w.velocity_m_s().powi(2))
            .sum();
        close(
            physical_change(flow.flow().gas(), before, flow.flow().totals())
                + store_change
                + kinetic
                + exterior_work,
            0.0,
            8e-8,
        );
    }
}

#[test]
fn g1_shortened_two_speed_trial_replays_native_sources_on_moved_volumes() {
    let mechanism = mechanism(40.0, 500.0, 0.0);
    let mut flow = closed(&mechanism, 4, chemistry());
    let mut replay = closed(&mechanism, 4, chemistry());
    let mut operator = molecular(flow.flow());
    let mut replay_operator = molecular(replay.flow());
    let mut walls = [2.0, 3.0].map(|mass| Piston::new(mass, 0.0).unwrap());
    let report = flow
        .step_with_inertial_walls_and_dissipation(
            &mut operator,
            &mut walls,
            Time::new(0.2),
            ReactiveStepControl::default(),
            DissipationControl::default(),
            wall_loads(1.0),
            wall_controls(),
            || false,
        )
        .unwrap();
    assert!(report.gas.coupling_retries > 0);
    let dt = report.gas.elapsed.value();
    let speeds = report.midpoint_velocity_m_s;
    replay_operator
        .admit_motion(&replay.accepted, &mut || false)
        .unwrap();
    let first = replay_operator
        .advance_moving_trial(
            &mut replay.accepted,
            Time::new(0.5 * dt),
            DissipationControl::default(),
            speeds,
            &mut || false,
        )
        .unwrap();
    let mut chemical = ChemistryReport::default();
    react_cells(
        &mut replay.accepted,
        &mut replay.reactor,
        &mut replay.amounts,
        0.5 * dt,
        &mut chemical,
        &mut || false,
    )
    .unwrap();
    let volumes = replay.flow().mesh().volumes().to_vec();
    let pressure = replay
        .accepted
        .step_with_moving_walls(
            Time::new(dt),
            ReactiveStepControl::default().transport,
            speeds,
            || false,
        )
        .unwrap();
    assert_eq!(pressure.elapsed.value(), dt);
    assert_ne!(replay.flow().mesh().volumes(), volumes.as_slice());
    react_cells(
        &mut replay.accepted,
        &mut replay.reactor,
        &mut replay.amounts,
        0.5 * dt,
        &mut chemical,
        &mut || false,
    )
    .unwrap();
    let second = replay_operator
        .advance_moving_trial(
            &mut replay.accepted,
            Time::new(0.5 * dt),
            DissipationControl::default(),
            speeds,
            &mut || false,
        )
        .unwrap();
    assert_eq!(flow.flow().mesh(), replay.flow().mesh());
    assert_eq!(flow.flow().time(), replay.flow().time());
    for (&actual, &expected) in flow.accepted.states.iter().zip(&replay.accepted.states) {
        close(actual, expected, 3e-12);
    }
    assert_ne!(first.before, second.before);
    let halves = report.gas.dissipation.unwrap();
    for (actual, expected) in [halves.first, halves.second]
        .into_iter()
        .zip([first, second])
    {
        close(
            actual.left_boundary_change.energy,
            expected.left_boundary_change.energy,
            3e-13,
        );
        close(
            actual.right_boundary_change.energy,
            expected.right_boundary_change.energy,
            3e-13,
        );
    }
}

fn thermal_fixture() -> (
    ReactingEuler1d,
    DissipationOperator<ConstantTransport>,
    PistonThermalEnds<'static>,
    [Piston; 2],
) {
    let flow = closed(&mechanism(8.0, 30.0, 0.0), 4, chemistry());
    let operator = molecular(flow.flow());
    let thermal = PistonThermalEnds::new(
        flow.flow(),
        [finite(3.0, 0.4, 1.0), finite(0.6, 0.7, 1.5)],
        [None, None],
        || false,
    )
    .unwrap();
    let walls = [0.2, 0.5].map(|mass| Piston::new(mass, 0.0).unwrap());
    (flow, operator, thermal, walls)
}

#[test]
fn g4_cancelled_molecular_heat_trial_preserves_both_walls_stores_and_moved_bindings() {
    let (mut reference, mut dr, mut tr, mut wr) = thermal_fixture();
    heat_step(&mut reference, &mut dr, &mut tr, &mut wr, || false).unwrap();
    let states = reference.accepted.states.clone();
    let mesh = reference.flow().mesh().clone();
    let time = reference.flow().time();
    let original_walls = wr;
    let contacts = *tr.contacts();
    let last_heat = tr.last_report();
    let mut polls = 0;
    let expected = heat_step(&mut reference, &mut dr, &mut tr, &mut wr, || {
        polls += 1;
        false
    })
    .unwrap();
    assert!(polls > 50);
    for cutoff in [1, polls / 2, polls] {
        let (mut flow, mut operator, mut thermal, mut walls) = thermal_fixture();
        heat_step(&mut flow, &mut operator, &mut thermal, &mut walls, || false).unwrap();
        let mut seen = 0;
        let error = heat_step(&mut flow, &mut operator, &mut thermal, &mut walls, || {
            seen += 1;
            seen == cutoff
        })
        .unwrap_err();
        assert_eq!(error, ReactiveError::Cancelled);
        assert_eq!(seen, cutoff);
        assert_eq!(flow.accepted.states, states);
        assert_eq!(flow.flow().mesh(), &mesh);
        assert_eq!(operator.mesh(), &mesh);
        assert_eq!(thermal.mesh(), &mesh);
        assert_eq!(flow.flow().time(), time);
        assert_eq!(walls, original_walls);
        assert_eq!(*thermal.contacts(), contacts);
        assert_eq!(thermal.last_report(), last_heat);
        let resumed =
            heat_step(&mut flow, &mut operator, &mut thermal, &mut walls, || false).unwrap();
        assert_eq!(resumed.gas.reactive.elapsed, expected.gas.reactive.elapsed);
        for (a, b) in walls.into_iter().zip(wr) {
            close(a.velocity_m_s(), b.velocity_m_s(), 2e-10);
        }
        for (&a, &b) in flow.accepted.states.iter().zip(&reference.accepted.states) {
            close(a, b, 2e-10);
        }
        assert_eq!(operator.mesh(), flow.flow().mesh());
        assert_eq!(thermal.mesh(), flow.flow().mesh());
        for (a, b) in thermal.contacts().iter().zip(tr.contacts()) {
            close(temperature(*a), temperature(*b), 2e-10);
        }
    }
}
