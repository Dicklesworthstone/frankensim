//! G1/G3/G4: actual finite gases, thermal walls and shared original inventories.
use super::super::tests::{controls, energy, gas, interface, near};
use super::*;
use crate::thermochemical_gas::energy::parcel::transfer::species::{
    effusion::channel::CylindricalMolecularGeometry, network::SpeciesNetworkError,
};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

const M: f64 = 0.028;
const AREA: f64 = 1e-6;
const DT: f64 = 0.001;

fn channel(length: f64, datum: f64) -> DiffuseMolecularChannel {
    DiffuseMolecularChannel::try_new(
        interface(datum),
        CylindricalMolecularGeometry::try_new(
            1e-8,
            length,
            1e-4,
            "synthetic fully accommodating diffuse circular channel",
        )
        .unwrap(),
    )
    .unwrap()
}
fn link(channel: &DiffuseMolecularChannel, temperature: f64) -> DiffuseChannelNetworkLink<'_> {
    DiffuseChannelNetworkLink {
        id: 7,
        parcels: [0, 1],
        channel,
        area_m2: AREA,
        wall_temperature_k: temperature,
    }
}
fn endpoints() -> BTreeMap<u64, GasEnergySolve> {
    BTreeMap::from([(0, energy()), (1, energy())])
}
fn incident_moles(t: f64, p: f64) -> f64 {
    DT * AREA * p / (2.0 * core::f64::consts::PI * M * R * t).sqrt()
}

#[test]
fn hot_and_cold_walls_change_actual_gas_temperature_and_close_total_entropy() {
    let gases = BTreeMap::from([
        (0, gas(400.0, 100_000.0, 0.0)),
        (1, gas(400.0, 100_000.0, 0.0)),
    ]);
    let original = gases.clone();
    let law = channel(1e-8, 0.0);
    let n = gases[&0].mass_kg() / M;
    for wall_temperature in [250.0, 900.0] {
        let result = step_candidate(
            &gases,
            GasTransportPaths::default(),
            &[link(&law, wall_temperature)],
            DT,
            &endpoints(),
            controls(),
            || false,
        )
        .unwrap();
        let heat = incident_moles(400.0, 100_000.0)
            * (1.0 - law.geometry().direct_probability())
            * 2.0
            * R
            * (wall_temperature - 400.0);
        let expected = 400.0 + heat / (1.5 * R * n);
        for id in [0, 1] {
            near(result.parcels()[&id].mass_kg(), gases[&id].mass_kg(), 1e-18);
            near(
                result.parcels()[&id].state().state().temperature,
                expected,
                3e-8,
            );
        }
        near(result.step().heat_into_walls_j, -2.0 * heat, 1e-12);
        assert!(result.step().energy_residual_j.abs() <= controls().exchange.energy_tolerance_j);
        assert!(result.step().total_entropy_change_j_k > 0.0);
        assert_eq!(
            result.step().exchange.boundaries.len(),
            2,
            "both reflected original inventories are spent"
        );
    }
    assert_eq!(gases, original);
}

#[test]
fn finite_temperatures_match_independent_direct_returned_and_accommodated_beams() {
    let law = channel(2e-8, 0.0);
    let gases = BTreeMap::from([
        (0, gas(600.0, 150_000.0, 0.0)),
        (1, gas(300.0, 100_000.0, 0.0)),
    ]);
    let result = step_candidate(
        &gases,
        GasTransportPaths::default(),
        &[link(&law, 800.0)],
        DT,
        &endpoints(),
        controls(),
        || false,
    )
    .unwrap();
    let tau = law.geometry().transmission_probability();
    let direct = law.geometry().direct_probability();
    let incident = [
        incident_moles(600.0, 150_000.0),
        incident_moles(300.0, 100_000.0),
    ];
    let source_t = [600.0, 300.0];
    for side in 0..2 {
        let n = gases[&(side as u64)].mass_kg() / M;
        let next_n = n + tau * (incident[1 - side] - incident[side]);
        let delta_u = 2.0
            * R
            * (-incident[side] * source_t[side]
                + (1.0 - tau) * incident[side] * 800.0
                + (tau - direct) * incident[1 - side] * 800.0
                + direct * incident[1 - side] * source_t[1 - side]);
        let next_t = (1.5 * R * n * source_t[side] + delta_u) / (1.5 * R * next_n);
        near(
            result.parcels()[&(side as u64)].mass_kg(),
            next_n * M,
            1e-18,
        );
        near(
            result.parcels()[&(side as u64)].state().state().temperature,
            next_t,
            3e-8,
        );
    }
    assert!(result.step().total_entropy_change_j_k > 0.0);
}

#[test]
fn zero_length_is_the_original_ballistic_transaction_and_wall_is_inactive() {
    let law = channel(0.0, 0.0);
    let gases = BTreeMap::from([
        (0, gas(600.0, 150_000.0, 0.0)),
        (1, gas(300.0, 100_000.0, 0.0)),
    ]);
    let result = step_candidate(
        &gases,
        GasTransportPaths::default(),
        &[link(&law, f64::NAN)],
        DT,
        &endpoints(),
        controls(),
        || false,
    )
    .unwrap();
    let native_links = [EffusionNetworkLink {
        id: 7,
        parcels: [0, 1],
        interface: law.mouth(),
        area_m2: AREA,
    }];
    let native = transport::step_candidate_with_paths(
        &gases,
        GasTransportPaths {
            effusion: &native_links,
            ..Default::default()
        },
        DT,
        &endpoints(),
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(result.parcels(), native.parcels());
    assert_eq!(result.step().exchange, native.step().exchange);
    assert_eq!(result.step().heat_into_walls_j, 0.0);
}

#[test]
fn reflected_molecules_share_the_original_inventory_allowance() {
    let law = channel(100e-8, 0.0);
    let gases = BTreeMap::from([
        (0, gas(400.0, 100_000.0, 0.0)),
        (1, gas(400.0, 100_000.0, 0.0)),
    ]);
    let original = gases.clone();
    let fraction = incident_moles(400.0, 100_000.0) * M / gases[&0].mass_kg();
    let mut limited = controls();
    limited.exchange.maximum_species_fraction = fraction * 0.5;
    assert!(law.geometry().transmission_probability() < 0.5);
    assert!(matches!(
        step_candidate(
            &gases,
            GasTransportPaths::default(),
            &[link(&law, 400.0)],
            DT,
            &endpoints(),
            limited,
            || false
        ),
        Err(GasTransportError::Exchange(
            SpeciesNetworkError::Withdrawal { .. }
        ))
    ));
    assert_eq!(gases, original);
}

#[test]
fn native_paths_and_channels_reconstruct_once_with_reference_independent_wall_heat() {
    let mut baseline = None;
    for datum in [0.0, 1e24, -1e24] {
        let law = channel(1e-8, datum);
        let gases = BTreeMap::from([
            (0, gas(600.0, 150_000.0, datum)),
            (1, gas(300.0, 100_000.0, datum)),
        ]);
        let native_links = [EffusionNetworkLink {
            id: 21,
            parcels: [0, 1],
            interface: law.mouth(),
            area_m2: AREA,
        }];
        let result = step_candidate(
            &gases,
            GasTransportPaths {
                effusion: &native_links,
                ..Default::default()
            },
            &[link(&law, 800.0)],
            DT,
            &endpoints(),
            controls(),
            || false,
        )
        .unwrap();
        assert_eq!(result.step().exchange.parcels.len(), 2);
        assert_eq!(result.step().transport.effusion.len(), 1);
        assert_eq!(result.step().channels.len(), 1);
        let state = [
            result.parcels()[&0].state().state().temperature,
            result.parcels()[&1].state().state().temperature,
            result.step().heat_into_walls_j,
        ];
        if let Some(expected) = baseline {
            for (actual, expected) in state.into_iter().zip(expected) {
                near(actual, expected, 3e-8);
            }
        } else {
            baseline = Some(state);
        }
        assert!(result.step().energy_residual_j.abs() <= controls().exchange.energy_tolerance_j);
    }
}

#[test]
fn closure_source_refusal_and_late_cancellation_leave_both_original_gases_unchanged() {
    let law = channel(1e-8, 0.0);
    let gases = BTreeMap::from([
        (0, gas(600.0, 150_000.0, 0.0)),
        (1, gas(300.0, 100_000.0, 0.0)),
    ]);
    let original = gases.clone();
    let closed = DiffuseChannelNetworkLink {
        area_m2: 0.0,
        ..link(&law, f64::NAN)
    };
    let result = step_candidate(
        &gases,
        GasTransportPaths::default(),
        &[closed],
        DT,
        &endpoints(),
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(result.parcels(), &gases);
    assert!(
        step_candidate(
            &gases,
            GasTransportPaths::default(),
            &[link(&law, 5000.0)],
            DT,
            &endpoints(),
            controls(),
            || false
        )
        .is_err()
    );
    let mut polls = 0usize;
    step_candidate(
        &gases,
        GasTransportPaths::default(),
        &[link(&law, 800.0)],
        DT,
        &endpoints(),
        controls(),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let mut count = 0usize;
    assert!(matches!(
        step_candidate(
            &gases,
            GasTransportPaths::default(),
            &[link(&law, 800.0)],
            DT,
            &endpoints(),
            controls(),
            || {
                count += 1;
                count == polls
            }
        ),
        Err(GasTransportError::Cancelled)
    ));
    assert_eq!(gases, original);
}
