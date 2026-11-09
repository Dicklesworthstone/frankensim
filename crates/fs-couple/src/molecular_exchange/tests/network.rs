//! G0/G1/G3/G4: one snapshot, physical moments, original shared inventory.
mod thermal;
mod trajectory;

use super::*;
use crate::molecular_exchange::network::{
    self, MolecularNetworkControls, MolecularNetworkError, MolecularNetworkLink,
};
use std::collections::BTreeMap;

fn budgets(frame: [f64; 3]) -> MolecularNetworkControls {
    let c = controls();
    MolecularNetworkControls {
        energy_frame_velocity_m_s: frame,
        maximum_duration_s: c.maximum_duration_s,
        maximum_species_fraction: c.maximum_species_fraction,
        mass_absolute_kg: c.mass_absolute_kg,
        momentum_absolute_n_s: c.momentum_absolute_n_s,
        energy_absolute_j: c.energy_absolute_j,
        relative_tolerance: c.relative_tolerance,
    }
}
fn solves(gases: &BTreeMap<u64, MovingGasParcel>) -> BTreeMap<u64, MovingGasBoundarySolve> {
    gases.keys().map(|&id| (id, controls().gas[0])).collect()
}
fn link<'a>(id: u64, gases: [u64; 2], law: &'a EffusiveDriftInterface) -> MolecularNetworkLink<'a> {
    MolecularNetworkLink {
        id,
        gases,
        interface: law,
        surface: surface(),
    }
}
fn three(model: &ThermochemicalGasModel, boost: [f64; 3]) -> BTreeMap<u64, MovingGasParcel> {
    [
        (2, 400.0, 180_000.0, [10.0, 30.0, -5.0]),
        (4, 700.0, 120_000.0, [-15.0, -20.0, 4.0]),
        (8, 500.0, 150_000.0, [5.0, 10.0, 2.0]),
    ]
    .into_iter()
    .map(|(id, t, p, v)| {
        (
            id,
            gas(
                model.clone(),
                t,
                p,
                core::array::from_fn(|d| v[d] + boost[d]),
            ),
        )
    })
    .collect()
}

#[test]
fn single_aperture_reproduces_the_original_pair_exactly() {
    let model = catalog(MIXTURE, 0.0);
    let law = interface(&model, 2000.0);
    let pair = [
        gas(model.clone(), 400.0, 180_000.0, [0.0, 30.0, 0.0]),
        gas(model, 700.0, 120_000.0, [0.0, -20.0, 0.0]),
    ];
    let expected = advance(&pair, &law, surface());
    let gases = BTreeMap::from([(2, pair[0].clone()), (4, pair[1].clone())]);
    let result = network::step_candidate(
        &gases,
        &[link(10, [2, 4], &law)],
        DT,
        &solves(&gases),
        budgets([0.0; 3]),
        || false,
    )
    .unwrap();
    assert_eq!(&result.gases()[&2], &expected.gases()[0]);
    assert_eq!(&result.gases()[&4], &expected.gases()[1]);
    assert_eq!(result.step().links[0].exchange, expected.step().exchange);
}

#[test]
fn three_gases_share_one_snapshot_and_conserve_actual_mass_momentum_and_energy() {
    let model = catalog(MIXTURE, 0.0);
    let law = interface(&model, 2000.0);
    let gases = three(&model, [0.0; 3]);
    let before = gases.clone();
    let mut links = [
        link(12, [2, 4], &law),
        link(7, [4, 8], &law),
        link(22, [8, 2], &law),
    ];
    links[1].surface.normal_left_to_right = [0.0, 1.0, 0.0];
    links[2].surface.normal_left_to_right = [0.0, 0.0, 1.0];
    links[0].surface.velocity_world_m_s = [2.0, -4.0, 1.0];
    links[1].surface.velocity_world_m_s = [-3.0, 5.0, 2.0];
    let result = network::step_candidate(
        &gases,
        &links,
        DT,
        &solves(&gases),
        budgets([7.0, -2.0, 4.0]),
        || false,
    )
    .unwrap();
    assert_eq!(gases, before);
    near(
        result.gases().values().map(|g| g.parcel().mass_kg()).sum(),
        gases.values().map(|g| g.parcel().mass_kg()).sum(),
        1e-15,
    );
    for d in 0..3 {
        near(
            result
                .gases()
                .values()
                .map(|g| g.momentum_kg_m_s()[d])
                .sum(),
            gases.values().map(|g| g.momentum_kg_m_s()[d]).sum(),
            1e-14,
        );
    }
    near(
        result
            .gases()
            .values()
            .map(|g| internal(g) + g.kinetic_energy_j())
            .sum(),
        gases
            .values()
            .map(|g| internal(g) + g.kinetic_energy_j())
            .sum(),
        2e-7,
    );
    assert!(result.step().pooled_stream_energy_residual_j.abs() < 1e-11);
    links.reverse();
    let reordered = network::step_candidate(
        &gases,
        &links,
        DT,
        &solves(&gases),
        budgets([7.0, -2.0, 4.0]),
        || false,
    )
    .unwrap();
    assert_eq!(result, reordered);
}

#[test]
fn shared_gross_withdrawal_cannot_be_funded_by_simultaneous_counterflow() {
    let model = catalog(PURE, 0.0);
    let law = interface(&model, 2000.0);
    let original = gas(model, T, P, [0.0; 3]);
    let gases = BTreeMap::from([
        (2, original.clone()),
        (4, original.clone()),
        (8, original.clone()),
    ]);
    let links = [link(1, [2, 4], &law), link(2, [4, 8], &law)];
    let mut c = budgets([0.0; 3]);
    c.maximum_species_fraction = 1.5 * crossing_mass(T, P, M) / original.parcel().mass_kg();
    for row in links {
        network::step_candidate(&gases, &[row], DT, &solves(&gases), c, || false).unwrap();
    }
    assert!(matches!(
        network::step_candidate(&gases, &links, DT, &solves(&gases), c, || false),
        Err(MolecularNetworkError::Withdrawal { gas: 4, .. })
    ));
}

#[test]
fn reversal_is_physical_and_multiple_apertures_do_not_create_sequential_ordering() {
    let model = catalog(PURE, 0.0);
    let law = interface(&model, 2000.0);
    let original = gas(model, T, P, [0.0; 3]);
    let gases = BTreeMap::from([(2, original.clone()), (4, original.clone()), (8, original)]);
    let mut links = [link(1, [2, 4], &law), link(2, [4, 8], &law)];
    let result = network::step_candidate(
        &gases,
        &links,
        DT,
        &solves(&gases),
        budgets([0.0; 3]),
        || false,
    )
    .unwrap();
    // The middle gas has exactly opposed equal crossing impulses, so one
    // simultaneous update leaves it at rest with its original thermal state.
    assert_eq!(&result.gases()[&4], &gases[&4]);
    for row in &mut links {
        row.gases.reverse();
        row.surface.normal_left_to_right = row.surface.normal_left_to_right.map(|v| -v);
    }
    let reversed = network::step_candidate(
        &gases,
        &links,
        DT,
        &solves(&gases),
        budgets([0.0; 3]),
        || false,
    )
    .unwrap();
    assert_eq!(result.gases(), reversed.gases());
}

#[test]
fn distinct_aperture_frames_and_large_formation_references_preserve_boosted_thermal_physics() {
    let evolve = |boost: [f64; 3], datum: f64| {
        let model = catalog(MIXTURE, datum);
        let law = interface(&model, 2000.0);
        let gases = three(&model, boost);
        let mut links = [link(1, [2, 4], &law), link(2, [4, 8], &law)];
        links[0].surface.velocity_world_m_s =
            core::array::from_fn(|d| [2.0, -4.0, 1.0][d] + boost[d]);
        links[1].surface.velocity_world_m_s =
            core::array::from_fn(|d| [-3.0, 5.0, 2.0][d] + boost[d]);
        links[1].surface.normal_left_to_right = [0.0, 1.0, 0.0];
        network::step_candidate(
            &gases,
            &links,
            DT,
            &solves(&gases),
            budgets(core::array::from_fn(|d| [7.0, -2.0, 4.0][d] + boost[d])),
            || false,
        )
        .unwrap()
    };
    let rest = evolve([0.0; 3], 0.0);
    let boost = [500.0, -700.0, 300.0];
    let moved = evolve(boost, 1e12);
    for (&id, gas) in rest.gases() {
        near(temperature(&moved.gases()[&id]), temperature(gas), 2e-8);
        near(
            moved.gases()[&id].parcel().mass_kg(),
            gas.parcel().mass_kg(),
            2e-15,
        );
        for d in 0..3 {
            near(
                moved.gases()[&id].velocity_m_s()[d] - boost[d],
                gas.velocity_m_s()[d],
                2e-10,
            );
        }
    }
}

#[test]
fn new_species_reach_all_original_receivers_without_forwarding_same_step_arrivals() {
    let model = catalog(MIXTURE, 0.0);
    let law = interface(&model, 2000.0);
    let only_a = model
        .with_composition(Composition::new(CompositionBasis::MoleFraction, vec![1.0, 0.0]).unwrap())
        .unwrap();
    let only_b = model
        .with_composition(Composition::new(CompositionBasis::MoleFraction, vec![0.0, 1.0]).unwrap())
        .unwrap();
    let gases = BTreeMap::from([
        (2, gas(only_a.clone(), T, P, [0.0; 3])),
        (4, gas(only_a, T, P, [0.0; 3])),
        (8, gas(only_b, T, P, [0.0; 3])),
    ]);
    let result = network::step_candidate(
        &gases,
        &[link(1, [2, 4], &law), link(2, [4, 8], &law)],
        DT,
        &solves(&gases),
        budgets([0.0; 3]),
        || false,
    )
    .unwrap();
    let b = SpeciesId::new("B").unwrap();
    assert!(
        !result.gases()[&2]
            .parcel()
            .species_inventory_kg()
            .unwrap()
            .contains_key(&b)
    );
    assert!(result.gases()[&4].parcel().species_inventory_kg().unwrap()[&b] > 0.0);
}

#[test]
fn topology_controls_idle_and_terminal_cancellation_are_atomic() {
    let model = catalog(PURE, 0.0);
    let law = interface(&model, 2000.0);
    let gases = three(&model, [0.0; 3]);
    let original = gases.clone();
    let endpoints = solves(&gases);
    let links = [link(1, [2, 4], &law), link(2, [4, 8], &law)];
    let idle =
        network::step_candidate(&gases, &[], DT, &endpoints, budgets([0.0; 3]), || false).unwrap();
    assert_eq!(idle.gases(), &gases);
    for bad in [
        vec![links[0], links[0]],
        vec![link(3, [2, 2], &law)],
        vec![link(3, [2, 99], &law)],
    ] {
        assert!(matches!(
            network::step_candidate(&gases, &bad, DT, &endpoints, budgets([0.0; 3]), || false),
            Err(MolecularNetworkError::Invalid(_))
        ));
    }
    let mut missing = endpoints.clone();
    missing.remove(&8);
    assert!(matches!(
        network::step_candidate(&gases, &links, DT, &missing, budgets([0.0; 3]), || false),
        Err(MolecularNetworkError::Invalid(_))
    ));
    let mut bad = budgets([0.0; 3]);
    bad.maximum_species_fraction = 1.0;
    assert!(matches!(
        network::step_candidate(&gases, &links, DT, &endpoints, bad, || false),
        Err(MolecularNetworkError::Invalid(_))
    ));
    let mut polls = 0;
    network::step_candidate(&gases, &links, DT, &endpoints, budgets([0.0; 3]), || {
        polls += 1;
        false
    })
    .unwrap();
    let mut remaining = polls;
    assert!(matches!(
        network::step_candidate(&gases, &links, DT, &endpoints, budgets([0.0; 3]), || {
            remaining -= 1;
            remaining == 0
        }),
        Err(MolecularNetworkError::Cancelled)
    ));
    assert_eq!(gases, original);
}
