//! G1/G3/G4: finite blowdown, thermal return, shared gross budgets and atomicity.
//! Synthetic NASA coefficients are equation fixtures, not empirical gas claims.
use super::super::super::tests::{controls, energy, gas, interface, near};
use super::super::super::{EffusionNetworkError, EffusionNetworkLink};
use super::super::tests::arrival_fixture;
use super::*;
use crate::thermochemical_gas::{
    ThermochemicalGasModel,
    energy::parcel::transfer::species::effusion::{
        EffusiveInterface, EffusiveSpeciesTransmission,
        channel::CylindricalMolecularGeometry,
        collision::{EffusiveCollisionModel, EffusiveCollisionPair},
    },
    tests::{point, species_with_energy_reference, transport_in_domain},
};
use fs_thermochem::{Composition, CompositionBasis, mixture::FrozenIdealGasMixtureModelV1};

const M: f64 = 0.028;
const AREA: f64 = 1e-6;
const DT: f64 = 0.001;

fn channel_from(mouth: EffusiveInterface, length_m: f64) -> DiffuseMolecularChannel {
    DiffuseMolecularChannel::try_new(
        mouth,
        CylindricalMolecularGeometry::try_new(
            1e-8,
            length_m,
            1e-4,
            "synthetic fully accommodating vacuum channel",
        )
        .unwrap(),
    )
    .unwrap()
}

fn channel(length_m: f64, datum: f64) -> DiffuseMolecularChannel {
    channel_from(interface(datum), length_m)
}

fn link(
    channel: &DiffuseMolecularChannel,
    wall_temperature_k: f64,
) -> VacuumChannelNetworkLink<'_> {
    VacuumChannelNetworkLink {
        id: 7,
        parcel: 0,
        channel,
        area_m2: AREA,
        wall_temperature_k,
    }
}

fn endpoints(gases: &BTreeMap<u64, ClosedGasParcel>) -> BTreeMap<u64, GasEnergySolve> {
    gases.keys().map(|&id| (id, energy())).collect()
}

fn incident_moles(t: f64, p: f64, m: f64) -> f64 {
    DT * AREA * p / (2.0 * core::f64::consts::PI * m * R * t).sqrt()
}

#[test]
fn monatomic_hot_cold_and_ambient_walls_match_three_independent_beam_energies() {
    let law = channel(2e-8, 0.0);
    let gases = BTreeMap::from([(0, gas(400.0, 100_000.0, 0.0))]);
    let original = gases.clone();
    let n0 = gases[&0].mass_kg() / M;
    let incident = incident_moles(400.0, 100_000.0, M);
    let tau = law.geometry().transmission_probability();
    let direct = law.geometry().direct_probability();
    for wall_t in [250.0, 400.0, 900.0] {
        let result = step_candidate(
            &gases,
            GasTransportPaths::default(),
            &[link(&law, wall_t)],
            DT,
            &endpoints(&gases),
            controls(),
            || false,
        )
        .unwrap();
        let returned_heat = (1.0 - tau) * incident * 2.0 * R * (wall_t - 400.0);
        let delta_u = -tau * incident * 2.0 * R * 400.0 + returned_heat;
        let next_n = n0 - tau * incident;
        let expected_t = (1.5 * R * n0 * 400.0 + delta_u) / (1.5 * R * next_n);
        let escape_energy = 2.0 * R * incident * (direct * 400.0 + (tau - direct) * wall_t);
        let wall_heat = -2.0 * R * incident * (1.0 - direct) * (wall_t - 400.0);
        near(result.parcels()[&0].mass_kg(), next_n * M, 1e-18);
        near(
            result.parcels()[&0].state().state().temperature,
            expected_t,
            3e-8,
        );
        near(
            result.step().audit.escaped_mass_kg,
            tau * incident * M,
            1e-18,
        );
        near(result.step().audit.escaped_energy_j, escape_energy, 1e-11);
        near(result.step().audit.heat_into_walls_j, wall_heat, 1e-11);
        near(
            result.step().outlets[0].wall.heat_into_gas_j,
            returned_heat,
            1e-11,
        );
        let spent = &result.step().exchange.boundaries[0].transfer;
        near(spent.outgoing_mass_kg, incident * M, 1e-18);
        near(spent.incoming_mass_kg, (1.0 - tau) * incident * M, 1e-18);
        near(spent.net_mass_into_gas_kg, -tau * incident * M, 1e-18);
        assert_eq!(result.step().exchange.parcels.len(), 1);
        assert!(
            result.step().audit.energy_residual_j.abs() <= controls().exchange.energy_tolerance_j
        );
    }
    assert_eq!(gases, original);
}

#[test]
fn zero_length_is_exactly_the_existing_true_vacuum_transaction() {
    let law = channel(0.0, 0.0);
    let gases = BTreeMap::from([(0, gas(600.0, 150_000.0, 0.0))]);
    let result = step_candidate(
        &gases,
        GasTransportPaths::default(),
        &[link(&law, f64::NAN)],
        DT,
        &endpoints(&gases),
        controls(),
        || false,
    )
    .unwrap();
    let native = thin_vacuum::step_candidate(
        &gases,
        GasTransportPaths::default(),
        &[VacuumEffusionLink {
            id: 7,
            parcel: 0,
            interface: law.mouth(),
            area_m2: AREA,
        }],
        DT,
        &endpoints(&gases),
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(result.parcels(), native.parcels());
    assert_eq!(result.step().exchange, native.step().exchange);
    assert_eq!(result.step().audit.heat_into_walls_j, 0.0);
    assert!(result.step().outlets[0].wall.reflected_species.is_empty());
}

#[test]
fn returning_molecules_spend_the_original_incident_allowance() {
    let law = channel(100e-8, 0.0);
    let gases = BTreeMap::from([(0, gas(400.0, 100_000.0, 0.0))]);
    let original = gases.clone();
    let fraction = incident_moles(400.0, 100_000.0, M) * M / gases[&0].mass_kg();
    let mut limited = controls();
    limited.exchange.maximum_species_fraction = 0.5 * fraction;
    assert!(law.geometry().transmission_probability() < 0.5);
    // The escaped fraction alone fits: it cannot stand in for all incident gas.
    assert!(
        thin_vacuum::step_candidate(
            &gases,
            GasTransportPaths::default(),
            &[VacuumEffusionLink {
                id: 7,
                parcel: 0,
                interface: law.mouth(),
                area_m2: law.transmitted_area_m2(AREA).unwrap(),
            }],
            DT,
            &endpoints(&gases),
            limited,
            || false,
        )
        .is_ok()
    );
    assert!(matches!(
        step_candidate(
            &gases,
            GasTransportPaths::default(),
            &[link(&law, 400.0)],
            DT,
            &endpoints(&gases),
            limited,
            || false,
        ),
        Err(GasTransportError::Exchange(
            SpeciesNetworkError::Withdrawal { .. }
        ))
    ));
    assert_eq!(gases, original);
}

#[test]
fn native_counterstreams_and_channel_returns_cannot_spend_separate_allowances() {
    let law = channel(2e-8, 0.0);
    let gases = BTreeMap::from([
        (0, gas(400.0, 100_000.0, 0.0)),
        (1, gas(400.0, 100_000.0, 0.0)),
    ]);
    let native = [EffusionNetworkLink {
        id: 99,
        parcels: [0, 1],
        interface: law.mouth(),
        area_m2: AREA,
    }];
    let fraction = incident_moles(400.0, 100_000.0, M) * M / gases[&0].mass_kg();
    let mut limited = controls();
    limited.exchange.maximum_species_fraction = 1.5 * fraction;
    assert!(
        step_candidate(
            &gases,
            GasTransportPaths::default(),
            &[link(&law, 400.0)],
            DT,
            &endpoints(&gases),
            limited,
            || false,
        )
        .is_ok()
    );
    assert!(
        transport::step_candidate_with_paths(
            &gases,
            GasTransportPaths {
                effusion: &native,
                ..Default::default()
            },
            DT,
            &endpoints(&gases),
            limited,
            || false,
        )
        .is_ok()
    );
    assert!(matches!(
        step_candidate(
            &gases,
            GasTransportPaths {
                effusion: &native,
                ..Default::default()
            },
            &[link(&law, 400.0)],
            DT,
            &endpoints(&gases),
            limited,
            || false,
        ),
        Err(GasTransportError::Exchange(
            SpeciesNetworkError::Withdrawal { .. }
        ))
    ));
}

#[test]
fn simultaneous_native_transport_and_wall_cooling_ignore_large_formation_datums() {
    let mut baseline = None;
    for datum in [0.0, 1e24, -1e24] {
        let law = channel(1e-8, datum);
        let gases = BTreeMap::from([
            (0, gas(600.0, 150_000.0, datum)),
            (1, gas(300.0, 100_000.0, datum)),
        ]);
        let native = [EffusionNetworkLink {
            id: 99,
            parcels: [0, 1],
            interface: law.mouth(),
            area_m2: AREA,
        }];
        let result = step_candidate(
            &gases,
            GasTransportPaths {
                effusion: &native,
                ..Default::default()
            },
            &[link(&law, 800.0)],
            DT,
            &endpoints(&gases),
            controls(),
            || false,
        )
        .unwrap();
        assert_eq!(result.step().exchange.parcels.len(), 2);
        assert_eq!(result.step().transport.effusion.len(), 1);
        assert_eq!(result.step().outlets.len(), 1);
        let values = [
            result.parcels()[&0].state().state().temperature,
            result.parcels()[&1].state().state().temperature,
            result.step().audit.heat_into_walls_j,
            result.step().audit.energy_residual_j,
        ];
        if let Some(expected) = baseline {
            for (actual, expected) in values.into_iter().zip(expected) {
                near(actual, expected, 3e-8);
            }
        } else {
            baseline = Some(values);
        }
        assert!(
            result.step().audit.energy_residual_j.abs() <= controls().exchange.energy_tolerance_j
        );
    }
}

#[test]
fn a_complete_window_audits_actual_storage_with_one_fixed_budget() {
    let law = channel(1e-8, 1e24);
    let before = BTreeMap::from([(0, gas(400.0, 100_000.0, 1e24))]);
    let first = step_candidate(
        &before,
        GasTransportPaths::default(),
        &[link(&law, 800.0)],
        DT,
        &endpoints(&before),
        controls(),
        || false,
    )
    .unwrap();
    let second = step_candidate(
        first.parcels(),
        GasTransportPaths::default(),
        &[link(&law, 800.0)],
        DT,
        &endpoints(&before),
        controls(),
        || false,
    )
    .unwrap();
    let mut outlets = first.step().outlets.clone();
    outlets.extend(second.step().outlets.clone());
    let audited = audit_window(&before, second.parcels(), &outlets, controls(), || false).unwrap();
    near(
        audited.escaped_mass_kg,
        first.step().audit.escaped_mass_kg + second.step().audit.escaped_mass_kg,
        1e-18,
    );
    assert!(audited.energy_residual_j.abs() <= controls().exchange.energy_tolerance_j);
    assert!(
        audit_window(&before, second.parcels(), &outlets[1..], controls(), || {
            false
        },)
        .is_err(),
        "an omitted earlier escape cannot be hidden by per-step balances"
    );
    let mut corrupted = outlets.clone();
    corrupted[0].wall.heat_into_wall_j = 0.0;
    assert!(
        audit_window(&before, second.parcels(), &corrupted, controls(), || false,).is_err(),
        "thermostat energy belongs to the same complete-window ledger"
    );
}

#[test]
fn a_binary_gas_fractionates_and_uses_each_species_actual_heat_capacity() {
    let rows = [("A", 0.002, 2.5), ("B", 0.028, 3.5)];
    let composition = [0.4, 0.6];
    let model = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            rows.iter()
                .map(|&(id, m, cp)| species_with_energy_reference(id, m, cp, 0.0, 0.0))
                .collect(),
            Composition::new(CompositionBasis::MoleFraction, composition.to_vec()).unwrap(),
        )
        .unwrap(),
        rows.iter()
            .map(|&(id, _, _)| transport_in_domain(id, 1.7e-5, 0.025, [200.0, 2000.0], [1.0, 1e7]))
            .collect(),
        None,
    )
    .unwrap();
    let gas = ClosedGasParcel::new(model.clone(), point(400.0, 100_000.0), 0.001).unwrap();
    let fixture = interface(0.0);
    let mouth = EffusiveInterface::try_new(
        model,
        rows.iter()
            .map(|&(id, _, _)| {
                EffusiveSpeciesTransmission::try_new(
                    SpeciesId::new(id).unwrap(),
                    1.0,
                    "synthetic unit molecular mouth",
                )
                .unwrap()
            })
            .collect(),
        fixture.support(),
        fixture.regime().clone(),
    )
    .unwrap();
    let law = channel_from(mouth, 2e-8);
    let gases = BTreeMap::from([(0, gas)]);
    let result = step_candidate(
        &gases,
        GasTransportPaths::default(),
        &[link(&law, 800.0)],
        DT,
        &endpoints(&gases),
        controls(),
        || false,
    )
    .unwrap();
    let n0 = 100_000.0 * 0.001 / (R * 400.0);
    let tau = law.geometry().transmission_probability();
    let mut next_moles = [0.0; 2];
    let mut expected_mass = 0.0;
    let mut expected_energy = 0.0;
    let mut expected_capacity = 0.0;
    for (index, &(_, mass, cp)) in rows.iter().enumerate() {
        let old_n = n0 * composition[index];
        let incident = incident_moles(400.0, 100_000.0 * composition[index], mass);
        next_moles[index] = old_n - tau * incident;
        expected_mass += mass * next_moles[index];
        expected_capacity += (cp - 1.0) * R * next_moles[index];
        expected_energy += (cp - 1.0) * R * old_n * 400.0 - tau * incident * (cp - 0.5) * R * 400.0
            + (1.0 - tau) * incident * (cp - 0.5) * R * (800.0 - 400.0);
        near(
            result.step().outlets[0].escaped_species[index].mass_kg,
            tau * incident * mass,
            1e-18,
        );
    }
    let after = &result.parcels()[&0];
    near(after.mass_kg(), expected_mass, 1e-18);
    near(
        after.state().state().temperature,
        expected_energy / expected_capacity,
        3e-8,
    );
    let fractions = after
        .model()
        .thermodynamics()
        .mole_composition()
        .fractions();
    near(
        fractions[0],
        next_moles[0] / next_moles.iter().sum::<f64>(),
        1e-14,
    );
    assert!(
        fractions[0] < composition[0],
        "the lighter species escapes preferentially"
    );
}

#[test]
fn zero_area_is_inactive_even_when_wall_and_channel_sources_are_inapplicable() {
    let law = channel(1e-8, 1e24);
    let gases = BTreeMap::from([(0, gas(400.0, 100_000.0, 0.0))]);
    let closed = VacuumChannelNetworkLink {
        area_m2: 0.0,
        ..link(&law, f64::NAN)
    };
    let result = step_candidate(
        &gases,
        GasTransportPaths::default(),
        &[closed],
        DT,
        &endpoints(&gases),
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(result.parcels(), &gases);
    assert!(result.step().exchange.boundaries.is_empty());
    assert!(result.step().outlets[0].escaped_species.is_empty());
    assert_eq!(result.step().audit.escaped_mass_kg, 0.0);
    assert_eq!(result.step().audit.heat_into_walls_j, 0.0);
}

#[test]
fn outlet_order_is_canonical_and_duplicate_ids_refuse() {
    let law = channel(1e-8, 0.0);
    let gases = BTreeMap::from([(0, gas(400.0, 100_000.0, 0.0))]);
    let first = link(&law, 800.0);
    let second = VacuumChannelNetworkLink {
        id: 99,
        wall_temperature_k: 250.0,
        ..first
    };
    let solve = |links: &[VacuumChannelNetworkLink<'_>]| {
        step_candidate(
            &gases,
            GasTransportPaths::default(),
            links,
            DT,
            &endpoints(&gases),
            controls(),
            || false,
        )
    };
    assert_eq!(
        solve(&[first, second]).unwrap(),
        solve(&[second, first]).unwrap()
    );
    assert!(solve(&[first, first]).is_err());
}

#[test]
fn source_domain_refusal_and_late_cancellation_never_mutate_original_gases() {
    let law = channel(1e-8, 0.0);
    let mismatched = channel(1e-8, 1e24);
    let gases = BTreeMap::from([(0, gas(400.0, 100_000.0, 0.0))]);
    let original = gases.clone();
    for outlet in [link(&law, 5000.0), link(&mismatched, 800.0)] {
        assert!(
            step_candidate(
                &gases,
                GasTransportPaths::default(),
                &[outlet],
                DT,
                &endpoints(&gases),
                controls(),
                || false,
            )
            .is_err()
        );
    }
    let mut polls = 0usize;
    let result = step_candidate(
        &gases,
        GasTransportPaths::default(),
        &[link(&law, 800.0)],
        DT,
        &endpoints(&gases),
        controls(),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let rebound = BTreeMap::from([(0, gas(400.0, 100_000.0, 1e24))]);
    assert!(
        result
            .step()
            .admit_endpoint_driving_forces(&rebound, || false)
            .is_err()
    );
    let mut visited = 0usize;
    assert!(matches!(
        step_candidate(
            &gases,
            GasTransportPaths::default(),
            &[link(&law, 800.0)],
            DT,
            &endpoints(&gases),
            controls(),
            || {
                visited += 1;
                visited == polls
            },
        ),
        Err(GasTransportError::Cancelled)
    ));
    assert_eq!(gases, original);
}

#[test]
fn endpoint_collision_check_uses_the_entire_channel_length() {
    let mouth = interface(0.0);
    let species = SpeciesId::new("A").unwrap();
    let sigma = 1.380_649e-23 * 400.0 / (100_000.0 * 2.0_f64.sqrt() * 5e-5);
    let collision = EffusiveCollisionModel::try_new(
        mouth.gas_catalog().clone(),
        vec![
            EffusiveCollisionPair::try_new(
                [species.clone(), species],
                sigma,
                "synthetic hard-sphere total cross section",
            )
            .unwrap(),
        ],
        mouth.support(),
        "synthetic collision support for finite channel endpoint",
    )
    .unwrap();
    let law = channel_from(mouth.with_collision_model(collision).unwrap(), 100e-8);
    let gases = BTreeMap::from([(0, gas(400.0, 100_000.0, 0.0))]);
    let result = step_candidate(
        &gases,
        GasTransportPaths::default(),
        &[link(&law, 400.0)],
        DT,
        &endpoints(&gases),
        controls(),
        || false,
    )
    .unwrap();
    let compressed = BTreeMap::from([(0, gas(400.0, 1_000_000.0, 0.0))]);
    // The thin mouth remains collisionless while the full bore no longer does.
    assert!(law.mouth().admit_parcel(&compressed[&0], || false).is_ok());
    assert!(
        result
            .step()
            .admit_endpoint_driving_forces(&compressed, || false)
            .is_err()
    );
}

#[test]
fn native_species_arrival_rechecks_active_wall_calorics_before_publication() {
    let b = SpeciesId::new("B").unwrap();
    for (b_max, length, area, wall_t, refuse) in [
        (600.0, 1e-8, AREA, 900.0, true),
        (1800.0, 1e-8, AREA, 900.0, false),
        (600.0, 1e-8, 0.0, f64::NAN, false),
        (600.0, 0.0, AREA, f64::NAN, false),
    ] {
        let (gases, law) = arrival_fixture(b_max, length);
        let original = gases.clone();
        let native = [EffusionNetworkLink {
            id: 99,
            parcels: [2, 0],
            interface: law.mouth(),
            area_m2: AREA,
        }];
        // Every gas solve stays inside BOTH source charts. A narrow-case
        // refusal must come from final B accommodation at the 900 K wall,
        // not from an unrelated caloric bracket outside B's 600 K support.
        let solves: BTreeMap<_, _> = gases
            .keys()
            .map(|&id| {
                let mut solve = energy();
                solve.maximum_temperature_k = 600.0;
                (id, solve)
            })
            .collect();
        let outlet = VacuumChannelNetworkLink {
            area_m2: area,
            ..link(&law, wall_t)
        };
        let result = step_candidate(
            &gases,
            GasTransportPaths {
                effusion: &native,
                ..Default::default()
            },
            &[outlet],
            DT,
            &solves,
            controls(),
            || false,
        );
        if refuse {
            assert!(
                matches!(
                    result,
                    Err(GasTransportError::Effusion(EffusionNetworkError::Face {
                        id: 7,
                        ..
                    }))
                ),
                "newly active B must support the retained wall temperature"
            );
        } else {
            let candidate = result.unwrap();
            assert!(
                candidate.parcels()[&0].species_mass_kg(&b).unwrap() > 0.0,
                "native transport must actually supply B to the outlet gas"
            );
            if area == 0.0 || length == 0.0 {
                assert_eq!(candidate.step().audit.heat_into_walls_j, 0.0);
            }
        }
        assert_eq!(gases, original);
    }
}
