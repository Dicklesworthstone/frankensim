//! Independent staged calorimetry, latent storage and whole-transaction replay.
use super::*;
use crate::multicomponent::dissipation::ConstantTransport;
use crate::multicomponent::network::dissipating::{
    DissipatingNetwork1d, DissipatingNetworkError, DissipatingStepControl,
};
use crate::multicomponent::network::dissipation::tests::law;
use crate::multicomponent::network::wall_exchange::SharedWallLink;

fn port(cell: usize, conductance_w_per_k: f64) -> SharedWallPort {
    SharedWallPort {
        gas: SharedWallGas::Branch { branch: 0, cell },
        conductance_w_per_k,
    }
}
fn link(a: usize, b: usize, conductance_w_per_k: f64) -> SharedWallLink {
    SharedWallLink {
        walls: [a, b],
        conductance_w_per_k,
    }
}
fn constant(t: f64, capacity: f64, port: SharedWallPort) -> SharedWall<'static> {
    SharedWall {
        storage: SharedWallStorage {
            thermal: WallThermal::Lumped {
                temperature_k: t,
                heat_capacity_j_per_k: capacity,
                conductance_w_per_k: 0.0,
            },
            heat_storage: None,
        },
        ports: vec![port],
    }
}
fn temperature(wall: &NetworkWallOperator<'_>, index: usize) -> f64 {
    let storage = &wall.shared_wall(index).unwrap().storage;
    if let Some(source) = &storage.heat_storage {
        source.temperature_k(storage.thermal).unwrap()
    } else {
        wall_temperature(WallContact {
            drag: WallDrag::None,
            thermal: storage.thermal,
        })
    }
}
fn pair(a: f64, b: f64, ca: f64, cb: f64, integrated: f64) -> [f64; 2] {
    let q = integrated * (a - b) / (1.0 + integrated / ca + integrated / cb);
    [a - q / ca, b + q / cb]
}

#[test]
fn gas_then_connected_wall_stages_match_independent_three_wall_calorimetry() {
    let g = gas(1, false, 0.0);
    let mut n = network(&g, &[0.0], true);
    let initial = fields(&n);
    let capacities = [0.7, 1.1, 0.9];
    let temperatures = [1.0, 4.0, 8.0];
    let films = [0.3, 0.7, 0.4];
    let mut wall = adiabatic(&n)
        .with_shared_walls(
            &n,
            vec![
                constant(temperatures[0], capacities[0], port(0, films[0])),
                constant(temperatures[1], capacities[1], port(1, films[1])),
                constant(
                    temperatures[2],
                    capacities[2],
                    SharedWallPort {
                        gas: SharedWallGas::Plenum(0),
                        conductance_w_per_k: films[2],
                    },
                ),
            ],
            SharedWallControl::default(),
            || false,
        )
        .unwrap()
        .with_shared_wall_links(vec![link(0, 1, 5.0), link(1, 2, 7.0)], || false)
        .unwrap();
    let dt = 0.2;
    let report = wall
        .advance_for(&mut n, Time::new(dt), controls(), || false)
        .unwrap();
    let after = fields(&n);
    let mut expected = [0.0; 3];
    let mut gas_heat = 0.0;
    for i in 0..3 {
        let response = dt * films[i] / (1.0 + dt * films[i] / 0.5);
        expected[i] =
            (capacities[i] * temperatures[i] + response * 2.0) / (capacities[i] + response);
        let q = response * (expected[i] - 2.0);
        gas_heat += q;
        let (old, new) = if i < 2 {
            (&initial.0[0][i], &after.0[0][i])
        } else {
            (&initial.1[0], &after.1[0])
        };
        assert_eq!(old.partial_densities, new.partial_densities);
        assert_eq!(old.momentum, new.momentum);
        close(
            g.primitive(new.as_cell(), || false).unwrap().temperature,
            2.0 + q / 0.5,
            2e-10,
        );
    }
    let gas_stage = expected;
    [expected[0], expected[1]] = pair(
        expected[0],
        expected[1],
        capacities[0],
        capacities[1],
        5.0 * dt,
    );
    [expected[1], expected[2]] = pair(
        expected[1],
        expected[2],
        capacities[1],
        capacities[2],
        7.0 * dt,
    );
    for (i, t) in expected.iter().enumerate() {
        close(temperature(&wall, i), *t, 2e-10);
    }
    let receipts = wall.shared_wall_link_reports().unwrap();
    assert_eq!(receipts.len(), 2);
    close(receipts[0].initial_temperature_k[0], gas_stage[0], 2e-10);
    assert_eq!(receipts[1].initial_state[0], receipts[0].endpoint_state[1]);
    for (link, receipt) in wall.shared_wall_links().unwrap().iter().zip(receipts) {
        assert_eq!(receipt.elapsed, Time::new(dt));
        close(
            receipt.heat_a_to_b_j,
            dt * link.conductance_w_per_k
                * (receipt.endpoint_temperature_k[0] - receipt.endpoint_temperature_k[1]),
            2e-11,
        );
        close(receipt.energy_residual_j, 0.0, 2e-11);
    }
    let wall_change: f64 = (0..3)
        .map(|i| capacities[i] * (temperature(&wall, i) - temperatures[i]))
        .sum();
    close(wall_change + gas_heat, 0.0, 2e-10);
    close(report.heat_to_gas_j, gas_heat, 2e-10);
    close(report.wall_energy_change_j, wall_change, 2e-10);
    assert_eq!(report.external_heat_to_gas_j, 0.0);
}

#[test]
fn connected_original_latent_inventory_exits_plateau_without_heating_gas_twice() {
    let g = gas(1, false, 0.0);
    let mut n = network(&g, &[0.0], false);
    let before = fields(&n);
    let latent =
        MaterialWallCaloric::from_thermal_law(Plateau, specimen(6.0), [0.0, 20.0]).unwrap();
    let mut wall = adiabatic(&n)
        .with_shared_walls(
            &n,
            vec![
                shared_storage(latent, vec![port(0, 0.0)]),
                constant(2.0, 3.0, port(1, 0.0)),
            ],
            SharedWallControl::default(),
            || false,
        )
        .unwrap()
        .with_shared_wall_links(vec![link(0, 1, 20.0)], || false)
        .unwrap();
    let first = wall
        .advance_for(&mut n, Time::new(0.05), controls(), || false)
        .unwrap();
    let q = 3.0 / (1.0 + 1.0 / 3.0);
    close(wall_energy(&wall), 6.0 - q, 2e-10);
    close(temperature(&wall, 0), 5.0, 1e-14);
    close(temperature(&wall, 1), 2.0 + q / 3.0, 2e-10);
    assert_eq!(first.heat_to_gas_j, 0.0);
    for _ in 0..16 {
        let report = wall
            .advance_for(&mut n, Time::new(0.1), controls(), || false)
            .unwrap();
        close(
            wall_energy(&wall) + 3.0 * temperature(&wall, 1),
            12.0,
            3e-10,
        );
        assert_eq!(report.heat_to_gas_j, 0.0);
        assert_eq!(report.external_heat_to_gas_j, 0.0);
    }
    assert!(wall_energy(&wall) < 3.0);
    close(temperature(&wall, 0), 3.5, 1e-5);
    close(temperature(&wall, 1), 3.5, 1e-5);
    assert_eq!(fields(&n), before);
    assert_eq!(
        wall.shared_wall(0)
            .unwrap()
            .storage
            .heat_storage
            .as_ref()
            .unwrap()
            .storage_basis(),
        ThermalStorageBasis::Enthalpy
    );
}

type Flow = DissipatingNetwork1d<ConstantTransport>;
fn composed(limit: usize) -> (Flow, NetworkWallOperator<'static>) {
    let g = gas(1, false, 0.0);
    let n = network(&g, &[0.0], false);
    let wall = adiabatic(&n)
        .with_shared_walls(
            &n,
            vec![
                constant(8.0, 0.7, port(0, 0.0)),
                constant(2.0, 1.1, port(1, 0.0)),
            ],
            SharedWallControl {
                maximum_source_evaluations: limit,
            },
            || false,
        )
        .unwrap()
        .with_shared_wall_links(vec![link(0, 1, 5.0)], || false)
        .unwrap();
    let flow = DissipatingNetwork1d::new(n, law(&g, 0.0, 0.0, 0.0, 0.0), vec![], || false).unwrap();
    (flow, wall)
}
fn flow_controls() -> DissipatingStepControl {
    DissipatingStepControl {
        maximum_coupling_retries: 0,
        ..DissipatingStepControl::default()
    }
}
fn states(flow: &Flow) -> Vec<Conserved> {
    (0..2)
        .map(|i| owned(flow.branches()[0].cell(i).unwrap()))
        .collect()
}

#[test]
fn both_link_halves_share_budget_and_preserve_prior_receipts_on_late_cancellation() {
    let (mut flow, mut wall) = composed(65_536);
    let first = flow
        .step_with_wall(
            &mut wall,
            Time::new(0.01),
            flow_controls(),
            controls(),
            || false,
        )
        .unwrap();
    let receipts = wall.shared_wall_link_reports().unwrap();
    assert_eq!(receipts.len(), 2);
    assert_eq!(receipts[0].elapsed.value(), first.elapsed.value() / 2.0);
    assert_eq!(receipts[1].initial_state, receipts[0].endpoint_state);
    let expected_half = pair(8.0, 2.0, 0.7, 1.1, 5.0 * first.elapsed.value() / 2.0);
    let expected_end = pair(
        expected_half[0],
        expected_half[1],
        0.7,
        1.1,
        5.0 * first.elapsed.value() / 2.0,
    );
    for i in 0..2 {
        close(temperature(&wall, i), expected_end[i], 2e-10);
    }
    let work = wall.last_shared_wall_report().unwrap().source_evaluations;
    let (mut exact, mut exact_wall) = composed(work);
    assert_eq!(
        exact
            .step_with_wall(
                &mut exact_wall,
                Time::new(0.01),
                flow_controls(),
                controls(),
                || false
            )
            .unwrap(),
        first
    );
    let (mut refused, mut refused_wall) = composed(work - 1);
    assert!(
        refused
            .step_with_wall(
                &mut refused_wall,
                Time::new(0.01),
                flow_controls(),
                controls(),
                || false
            )
            .is_err()
    );
    assert_eq!(refused.time(), Time::new(0.0));
    assert!(refused_wall.shared_wall_link_reports().is_none());
    close(temperature(&refused_wall, 0), 8.0, 1e-14);
    // Work may differ at the next physical state; cancellation uses a fresh
    // generous allowance rather than the previous interval's exact cost.
    let (mut exact, mut exact_wall) = composed(65_536);
    assert_eq!(
        exact
            .step_with_wall(
                &mut exact_wall,
                Time::new(0.01),
                flow_controls(),
                controls(),
                || false
            )
            .unwrap(),
        first
    );
    let prior_states = states(&exact);
    let prior_receipts = exact_wall.shared_wall_link_reports().unwrap().to_vec();
    let prior_report = exact_wall.last_shared_wall_report();
    // Use a reference second interval to locate the last complete publication poll.
    let mut polls = 0;
    let expected = flow
        .step_with_wall(
            &mut wall,
            Time::new(0.01),
            flow_controls(),
            controls(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut seen = 0;
    let error = exact
        .step_with_wall(
            &mut exact_wall,
            Time::new(0.01),
            flow_controls(),
            controls(),
            || {
                seen += 1;
                seen >= polls
            },
        )
        .unwrap_err();
    assert_eq!(error, DissipatingNetworkError::Cancelled);
    assert_eq!(states(&exact), prior_states);
    assert_eq!(exact.time(), first.elapsed);
    assert_eq!(exact.last_report(), Some(first));
    assert_eq!(
        exact_wall.shared_wall_link_reports().unwrap(),
        prior_receipts
    );
    assert_eq!(exact_wall.last_shared_wall_report(), prior_report);
    assert_eq!(
        exact
            .step_with_wall(
                &mut exact_wall,
                Time::new(0.01),
                flow_controls(),
                controls(),
                || false
            )
            .unwrap(),
        expected
    );
    assert_eq!(
        exact_wall.shared_wall_link_reports(),
        wall.shared_wall_link_reports()
    );
}

#[test]
fn refused_link_restores_prior_gas_heat_and_original_wall_domain() {
    let g = gas(1, false, 0.0);
    let mut n = network(&g, &[0.0], false);
    let source =
        MaterialWallCaloric::from_thermal_law(UnitInternalStorage, specimen(5.0), [4.99, 6.0])
            .unwrap();
    let mut wall = adiabatic(&n)
        .with_shared_walls(
            &n,
            vec![
                shared_storage(source, vec![port(0, 0.001)]),
                constant(1.0, 1.0, port(1, 0.2)),
            ],
            SharedWallControl::default(),
            || false,
        )
        .unwrap()
        .with_shared_wall_links(vec![link(0, 1, 100.0)], || false)
        .unwrap();
    let before = fields(&n);
    let original = [
        wall.shared_wall(0).unwrap().clone(),
        wall.shared_wall(1).unwrap().clone(),
    ];
    assert!(
        wall.advance_for(&mut n, Time::new(0.1), controls(), || false)
            .is_err()
    );
    assert_eq!(fields(&n), before);
    for (i, state) in original.iter().enumerate() {
        assert_eq!(wall.shared_wall(i).unwrap(), state);
    }
    assert!(wall.shared_wall_link_reports().is_none());
    assert!(wall.last_report().is_none());
}
