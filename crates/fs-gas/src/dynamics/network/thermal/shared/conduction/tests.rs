//! G0/G1: original inventories, physical port equations and time refinement.

use super::super::super::shared_tests::{
    attach, base, close, controls, fields, gas, latent, lumped, periodic, port, step,
    wall_temperature,
};
use super::*;
use crate::dynamics::StepControl;

fn link(a: usize, b: usize, g: f64) -> SharedWallLink {
    SharedWallLink {
        walls: [a, b],
        conductance_w_per_k: g,
    }
}

fn network() -> GasNetwork1d {
    let gas = gas(0.0, 0.0);
    GasNetwork1d::new(vec![periodic(&gas, 3.0, 0.0)], vec![], vec![], || false).unwrap()
}

fn walls(
    network: &GasNetwork1d,
    temperatures: &[f64],
    links: Vec<SharedWallLink>,
) -> NetworkWallOperator<'static> {
    attach(
        network,
        WallDrag::None,
        temperatures
            .iter()
            .enumerate()
            .map(|(i, &t)| lumped(t, 1.0, vec![port(0, i, 0.0)]))
            .collect(),
    )
    .with_shared_wall_links(links, || false)
    .unwrap()
}

fn profile(walls: &NetworkWallOperator<'_>, n: usize) -> Vec<f64> {
    (0..n).map(|i| wall_temperature(walls, i)).collect()
}

fn saved<'law>(walls: &NetworkWallOperator<'law>, n: usize) -> Vec<SharedWall<'law>> {
    (0..n)
        .map(|i| walls.shared_wall(i).unwrap().clone())
        .collect()
}

fn pair(ta: f64, tb: f64, ca: f64, cb: f64, integrated: f64) -> [f64; 2] {
    let heat = integrated * (ta - tb) / (1.0 + integrated / ca + integrated / cb);
    [ta - heat / ca, tb + heat / cb]
}

#[test]
fn gas_stage_and_connected_wall_stages_match_independent_calorimetry() {
    let mut network = network();
    let capacities = [0.7, 1.1, 0.9];
    let initial = [1.0, 4.0, 2.0];
    let films = [2.0, 1.0, 3.0];
    let mut walls = attach(
        &network,
        WallDrag::None,
        (0..3)
            .map(|i| lumped(initial[i], capacities[i], vec![port(0, i, films[i])]))
            .collect(),
    )
    .with_shared_wall_links(vec![link(0, 1, 5.0), link(1, 2, 7.0)], || false)
    .unwrap();
    let before = network.totals().energy
        + initial
            .iter()
            .zip(capacities)
            .map(|(t, c)| t * c)
            .sum::<f64>();
    let report = step(&mut network, &mut walls, 0.01);
    let dt = report.elapsed.value();
    let gas_capacity = 0.175 * 2.5;
    let mut expected = [0.0; 3];
    let mut heat_to_gas = 0.0;
    for i in 0..3 {
        let a = dt * films[i];
        let response = a * gas_capacity / (a + gas_capacity);
        expected[i] = (capacities[i] * initial[i] + response * 3.0) / (capacities[i] + response);
        let tg = (gas_capacity * 3.0 + a * expected[i]) / (gas_capacity + a);
        close(
            network
                .gas()
                .primitive(network.branches()[0].states()[i], || false)
                .unwrap()
                .temperature,
            tg,
            2e-10,
        );
        heat_to_gas += gas_capacity * (tg - 3.0);
    }
    let gas_stage = expected;
    let first = pair(
        expected[0],
        expected[1],
        capacities[0],
        capacities[1],
        5.0 * dt,
    );
    expected[0] = first[0];
    expected[1] = first[1];
    let second = pair(
        expected[1],
        expected[2],
        capacities[1],
        capacities[2],
        7.0 * dt,
    );
    expected[1] = second[0];
    expected[2] = second[1];
    for (a, b) in profile(&walls, 3).iter().zip(expected) {
        close(*a, b, 2e-10);
    }
    let links = walls.shared_wall_links().unwrap();
    let receipts = walls.shared_wall_link_reports().unwrap();
    assert_eq!(receipts[1].initial_state[0], receipts[0].endpoint_state[1]);
    close(receipts[0].initial_temperature_k[0], gas_stage[0], 2e-10);
    assert_ne!(
        gas_stage[0], expected[0],
        "gas and link stages must be distinguished"
    );
    for (port, receipt) in links.iter().zip(receipts) {
        assert_eq!(receipt.elapsed, report.elapsed);
        close(receipt.energy_residual_j, 0.0, 2e-11);
        close(
            receipt.heat_a_to_b_j,
            dt * port.conductance_w_per_k
                * (receipt.endpoint_temperature_k[0] - receipt.endpoint_temperature_k[1]),
            2e-11,
        );
    }
    let source = report.thermal.unwrap();
    close(source.source.heat_to_gas_j, heat_to_gas, 2e-10);
    close(source.source.wall_energy_change_j, -heat_to_gas, 2e-10);
    close(source.gas_wall_conservation_residual_j, 0.0, 2e-11);
    close(
        network.totals().energy
            + profile(&walls, 3)
                .iter()
                .zip(capacities)
                .map(|(t, c)| t * c)
                .sum::<f64>(),
        before,
        2e-10,
    );
}

#[test]
fn three_wall_chain_converges_to_continuous_diffusion_and_loses_order_bias() {
    let finish: f64 = 0.2;
    let slow = (-4.0 * finish).exp();
    let fast = (-12.0 * finish).exp();
    let exact = [
        7.0 / 3.0 - 1.5 * slow + fast / 6.0,
        7.0 / 3.0 - fast / 3.0,
        7.0 / 3.0 + 1.5 * slow + fast / 6.0,
    ];
    let mut errors = Vec::new();
    let mut biases = Vec::new();
    for steps in [20, 40, 80] {
        let mut outcomes = Vec::new();
        for reversed in [false, true] {
            let mut network = network();
            let mut links = vec![link(0, 1, 4.0), link(1, 2, 4.0)];
            if reversed {
                links.reverse();
            }
            let mut walls = walls(&network, &[1.0, 2.0, 4.0], links);
            let original_gas = fields(&network);
            let mut entropy = 8.0_f64.ln();
            for _ in 0..steps {
                let report = step(&mut network, &mut walls, finish / steps as f64);
                close(report.elapsed.value(), finish / steps as f64, 1e-15);
                let t = profile(&walls, 3);
                close(t.iter().sum(), 7.0, 2e-10);
                assert!(t.iter().all(|t| (1.0..=4.0).contains(t)));
                let next = t.iter().map(|t| t.ln()).sum::<f64>();
                assert!(next >= entropy - 1e-11);
                entropy = next;
            }
            assert_eq!(fields(&network), original_gas);
            outcomes.push(profile(&walls, 3));
        }
        errors.push(
            outcomes[0]
                .iter()
                .zip(exact)
                .map(|(a, b)| (a - b).abs())
                .fold(0.0_f64, f64::max),
        );
        biases.push(
            outcomes[0]
                .iter()
                .zip(&outcomes[1])
                .map(|(a, b)| (a - b).abs())
                .fold(0.0_f64, f64::max),
        );
    }
    for values in [errors, biases] {
        assert!(values[0] > 1e-4, "nontrivial composition error {values:?}");
        for pair in values.windows(2) {
            assert!(pair[1] < 0.6 * pair[0], "first-order refinement {values:?}");
        }
    }
}

#[test]
fn latent_inventory_transfers_heat_on_plateau_then_crosses_into_sensible_region() {
    let mut network = network();
    let mut walls = attach(
        &network,
        WallDrag::None,
        vec![
            lumped(1.0, 1.0, vec![port(0, 0, 0.0)]),
            latent([0.0, 10.0], vec![port(0, 1, 0.0)]),
        ],
    )
    .with_shared_wall_links(vec![link(1, 0, 100.0)], || false)
    .unwrap();
    let initial_source = walls.shared_wall(1).unwrap().storage.heat_storage.clone();
    let report = step(&mut network, &mut walls, 0.001);
    let a = report.elapsed.value() * 100.0;
    let q = a * 2.0 / (1.0 + a);
    close(wall_temperature(&walls, 0), 1.0 + q, 2e-10);
    assert_eq!(wall_temperature(&walls, 1), 3.0);
    let WallThermal::Caloric {
        specific_energy_j_kg,
        mass_kg,
        ..
    } = walls.shared_wall(1).unwrap().storage.thermal
    else {
        panic!("original h state")
    };
    assert_eq!(mass_kg, 0.2);
    close(specific_energy_j_kg, 4.0 - q / 0.2, 3e-10);
    assert!((2.0..5.0).contains(&specific_energy_j_kg));
    for _ in 0..16 {
        step(&mut network, &mut walls, 0.01);
        let WallThermal::Caloric {
            specific_energy_j_kg,
            mass_kg,
            ..
        } = walls.shared_wall(1).unwrap().storage.thermal
        else {
            panic!("original h state")
        };
        close(
            wall_temperature(&walls, 0) + mass_kg * specific_energy_j_kg,
            1.8,
            2e-10,
        );
        assert_eq!(
            walls.shared_wall(1).unwrap().storage.heat_storage,
            initial_source
        );
    }
    close(wall_temperature(&walls, 0), 5.0 / 3.0, 2e-8);
    close(wall_temperature(&walls, 1), 5.0 / 3.0, 2e-8);
}

#[test]
fn cyclic_parallel_ports_preserve_exact_equilibrium_and_isolated_storage() {
    let mut network = network();
    let mut walls = walls(
        &network,
        &[3.0, 3.0, 3.0, 4.0],
        vec![
            link(0, 1, 1e12),
            link(1, 2, 3e12),
            link(2, 0, 4e12),
            link(0, 1, 2e12),
            link(0, 3, 0.0),
        ],
    );
    let original = saved(&walls, 4);
    let original_gas = fields(&network);
    for _ in 0..4 {
        let report = step(&mut network, &mut walls, 0.01);
        assert_eq!(saved(&walls, 4), original);
        assert_eq!(fields(&network), original_gas);
        assert_eq!(report.thermal.unwrap().source.wall_energy_change_j, 0.0);
        for link in walls.shared_wall_link_reports().unwrap() {
            assert_eq!(link.initial_state, link.endpoint_state);
            assert_eq!(link.heat_a_to_b_j, 0.0);
            assert_eq!(link.energy_residual_j, 0.0);
            assert_eq!(link.flux_residual_j, 0.0);
        }
    }
}

#[test]
fn invalid_links_domain_and_whole_attempt_work_limits_refuse_atomically() {
    let mut network = network();
    let make = |network: &GasNetwork1d| {
        walls(
            network,
            &[1.0, 4.0, 2.0],
            vec![link(0, 1, 3.0), link(1, 2, 2.0)],
        )
    };
    for bad in [
        link(0, 0, 1.0),
        link(0, 3, 1.0),
        link(0, 1, -1.0),
        link(0, 1, f64::NAN),
    ] {
        assert!(
            attach(
                &network,
                WallDrag::None,
                vec![
                    lumped(1.0, 1.0, vec![port(0, 0, 0.0)]),
                    lumped(4.0, 1.0, vec![port(0, 1, 0.0)]),
                ]
            )
            .with_shared_wall_links(vec![bad], || false)
            .is_err()
        );
    }
    assert!(
        base(&network, WallDrag::None)
            .with_shared_wall_links(vec![link(0, 1, 1.0)], || false)
            .is_err()
    );
    let mut good = make(&network);
    step(&mut network, &mut good, 0.01);
    let needed = good.last_shared_wall_report().unwrap().source_evaluations;
    assert!(needed > 10);
    let mut exact_network = super::tests::network();
    let mut exact = make(&exact_network);
    exact
        .shared
        .as_mut()
        .unwrap()
        .controls
        .maximum_source_evaluations = needed;
    step(&mut exact_network, &mut exact, 0.01);
    assert_eq!(saved(&good, 3), saved(&exact, 3));
    let mut refused_network = super::tests::network();
    let mut refused = make(&refused_network);
    refused
        .shared
        .as_mut()
        .unwrap()
        .controls
        .maximum_source_evaluations = needed - 1;
    let original = saved(&refused, 3);
    let original_gas = fields(&refused_network);
    let no_retry = StepControl {
        maximum_halvings: 0,
        ..StepControl::default()
    };
    assert!(
        refused_network
            .step_with_walls(Time::new(0.01), no_retry, &mut refused, controls(), || {
                false
            })
            .is_err()
    );
    assert_eq!(saved(&refused, 3), original);
    assert_eq!(fields(&refused_network), original_gas);
    assert_eq!(refused.shared_wall_link_reports(), None);
    assert_eq!(refused_network.time(), Time::new(0.0));

    // First link succeeds in scratch; the last bounded source cannot accept heat.
    let mut domain = attach(
        &refused_network,
        WallDrag::None,
        vec![
            lumped(4.0, 1.0, vec![port(0, 0, 0.0)]),
            lumped(4.5, 1.0, vec![port(0, 1, 0.0)]),
            latent([4.0, 4.000001], vec![port(0, 2, 0.0)]),
        ],
    )
    .with_shared_wall_links(vec![link(0, 1, 2.0), link(1, 2, 3.0)], || false)
    .unwrap();
    let original = saved(&domain, 3);
    assert!(
        refused_network
            .step_with_walls(Time::new(0.01), no_retry, &mut domain, controls(), || false)
            .is_err()
    );
    assert_eq!(saved(&domain, 3), original);
    assert_eq!(fields(&refused_network), original_gas);
    assert_eq!(domain.shared_wall_link_reports(), None);
}

#[test]
fn final_cancellation_preserves_link_history_and_retries_bitwise() {
    let make = || {
        let network = network();
        let walls = walls(
            &network,
            &[1.0, 4.0, 2.0],
            vec![link(0, 1, 3.0), link(1, 2, 2.0)],
        );
        (network, walls)
    };
    let (mut good, mut good_walls) = make();
    let (mut cancelled, mut cancelled_walls) = make();
    step(&mut good, &mut good_walls, 0.01);
    step(&mut cancelled, &mut cancelled_walls, 0.01);
    let states = saved(&cancelled_walls, 3);
    let old_gas = fields(&cancelled);
    let old_time = cancelled.time();
    let old_report = cancelled.last_report();
    let old_shared = cancelled_walls.last_shared_wall_report();
    let old_links = cancelled_walls.shared_wall_link_reports().unwrap().to_vec();
    let mut calls = 0;
    good.step_with_walls(
        Time::new(0.01),
        StepControl::default(),
        &mut good_walls,
        controls(),
        || {
            calls += 1;
            false
        },
    )
    .unwrap();
    let mut count = 0;
    assert_eq!(
        cancelled.step_with_walls(
            Time::new(0.01),
            StepControl::default(),
            &mut cancelled_walls,
            controls(),
            || {
                count += 1;
                count == calls
            }
        ),
        Err(Error::Cancelled)
    );
    assert_eq!(saved(&cancelled_walls, 3), states);
    assert_eq!(fields(&cancelled), old_gas);
    assert_eq!(cancelled.time(), old_time);
    assert_eq!(cancelled.last_report(), old_report);
    assert_eq!(cancelled_walls.last_shared_wall_report(), old_shared);
    assert_eq!(
        cancelled_walls.shared_wall_link_reports().unwrap(),
        old_links
    );
    step(&mut cancelled, &mut cancelled_walls, 0.01);
    assert_eq!(saved(&cancelled_walls, 3), saved(&good_walls, 3));
    assert_eq!(fields(&cancelled), fields(&good));
    assert_eq!(cancelled.last_report(), good.last_report());
    assert_eq!(
        cancelled_walls.shared_wall_link_reports(),
        good_walls.shared_wall_link_reports()
    );
}
