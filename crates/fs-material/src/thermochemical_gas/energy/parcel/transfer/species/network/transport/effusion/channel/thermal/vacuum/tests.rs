//! G1/G3/G4: independent beam/storage balances and atomic finite-wall feedback.
use super::super::super::super::tests::{controls, gas, interface, near};
use super::super::tests::{channel, endpoints, material};
use super::*;
use crate::thermal::transport::volumetric::timeline::FixedVolumeThermalLaw;
use crate::thermochemical_gas::energy::parcel::transfer::species::{
    effusion::channel::{ClausingControls, CylindricalMolecularGeometry},
    network::{SpeciesNetworkError, transport::effusion::EffusionNetworkLink},
};
use fs_thermal::enthalpy::ThermalLaw;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

const AREA: f64 = 1e-6;
const DT: f64 = 0.001;
const M: f64 = 0.028;

fn link(channel: &DiffuseMolecularChannel) -> FiniteVacuumChannelLink<'_> {
    FiniteVacuumChannelLink {
        id: 7,
        parcel: 0,
        channel,
        area_m2: AREA,
        wall: 9,
    }
}

fn incident_moles(t: f64, p: f64) -> f64 {
    DT * AREA * p / (2.0 * core::f64::consts::PI * M * R * t).sqrt()
}

#[test]
fn partial_accommodation_changes_finite_wall_storage_and_escaped_energy_once() {
    let source = material(1e6, 0.0, [200.0, 1800.0]);
    let law = FixedVolumeThermalLaw::try_new(&source, 1000.0, 200.0)
        .unwrap()
        .with_sensible_entropy();
    let walls = BTreeMap::from([(
        9,
        FixedVolumeChannelWall::try_new_isochoric(&law, 1e-6, 900.0).unwrap(),
    )]);
    let path = channel(2e-8)
        .with_thermal_accommodation(
            0.25,
            ClausingControls {
                initial_rings: 8,
                maximum_rings: 512,
                maximum_iterations: 256,
                absolute_tolerance: 1e-6,
                relative_tolerance: 1e-4,
            },
            "synthetic shared per-collision energy accommodation",
            || false,
        )
        .unwrap();
    let gases = BTreeMap::from([(0, gas(400.0, 100_000.0, 0.0))]);
    let result = step_candidate(
        &gases,
        &walls,
        GasTransportPaths::default(),
        &[link(&path)],
        DT,
        &endpoints(&gases),
        controls(),
        || false,
    )
    .unwrap();
    let model = path.thermal_accommodation().unwrap();
    let incident = incident_moles(400.0, 100_000.0);
    let tau = path.geometry().transmission_probability();
    let returned_heat = incident * model.reflected_heat_probability() * 2.0 * R * 500.0;
    let escaped_heat = incident * model.transmitted_heat_probability() * 2.0 * R * 500.0;
    let n0 = gases[&0].mass_kg() / M;
    let n1 = n0 - tau * incident;
    let tg =
        (1.5 * R * n0 * 400.0 - tau * incident * 2.0 * R * 400.0 + returned_heat) / (1.5 * R * n1);
    near(result.parcels()[&0].state().state().temperature, tg, 3e-8);
    near(
        result.step().gas.audit.escaped_energy_j,
        tau * incident * 2.0 * R * 400.0 + escaped_heat,
        1e-11,
    );
    near(
        result.walls()[&9].temperature_k(),
        900.0 - returned_heat - escaped_heat,
        1e-10,
    );
    near(
        result.step().wall_energy_change_j,
        -returned_heat - escaped_heat,
        1e-9,
    );
    assert!(result.step().energy_residual_j.abs() <= controls().exchange.energy_tolerance_j);
}

#[test]
fn elastic_channel_keeps_finite_wall_unchanged_outside_gas_caloric_support() {
    let source = material(1e6, 0.0, [1800.0, 2000.0]);
    let law = FixedVolumeThermalLaw::try_new(&source, 1000.0, 1800.0)
        .unwrap()
        .with_sensible_entropy();
    let walls = BTreeMap::from([(
        9,
        FixedVolumeChannelWall::try_new_isochoric(&law, 1e-6, 1900.0).unwrap(),
    )]);
    let path = channel(2e-8)
        .with_thermal_accommodation(
            0.0,
            ClausingControls::default(),
            "synthetic exactly elastic wall",
            || false,
        )
        .unwrap();
    let gases = BTreeMap::from([(0, gas(400.0, 100_000.0, 0.0))]);
    let result = step_candidate(
        &gases,
        &walls,
        GasTransportPaths::default(),
        &[link(&path)],
        DT,
        &endpoints(&gases),
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(
        result.walls()[&9].specific_energy_j_kg(),
        walls[&9].specific_energy_j_kg()
    );
    assert_eq!(result.step().wall_energy_change_j, 0.0);
    assert_eq!(result.step().wall_entropy_change_j_k, 0.0);
    assert!(result.step().gas.outlets[0].wall.reflected_species[0].mass_kg > 0.0);
    near(
        result.step().gas.exchange.boundaries[0]
            .transfer
            .outgoing_mass_kg,
        incident_moles(400.0, 100_000.0) * M,
        1e-18,
    );
    assert!(result.parcels()[&0].mass_kg() < gases[&0].mass_kg());
}

#[test]
fn hot_and_cold_finite_walls_match_independent_gas_escape_and_nonlinear_storage() {
    let source = material(1e6, 1000.0, [200.0, 1800.0]);
    let law = FixedVolumeThermalLaw::try_new(&source, 1000.0, 200.0)
        .unwrap()
        .with_sensible_entropy();
    let path = channel(2e-8);
    let gases = BTreeMap::from([(0, gas(400.0, 100_000.0, 0.0))]);
    let n0 = gases[&0].mass_kg() / M;
    let incident = incident_moles(400.0, 100_000.0);
    let tau = path.geometry().transmission_probability();
    let direct = path.geometry().direct_probability();
    for tw in [250.0, 900.0] {
        let walls = BTreeMap::from([(
            9,
            FixedVolumeChannelWall::try_new_isochoric(&law, 1e-6, tw).unwrap(),
        )]);
        let result = step_candidate(
            &gases,
            &walls,
            GasTransportPaths::default(),
            &[link(&path)],
            DT,
            &endpoints(&gases),
            controls(),
            || false,
        )
        .unwrap();
        let returned_heat = (1.0 - tau) * incident * 2.0 * R * (tw - 400.0);
        let delta_u = -tau * incident * 2.0 * R * 400.0 + returned_heat;
        let n1 = n0 - tau * incident;
        let tg = (1.5 * R * n0 * 400.0 + delta_u) / (1.5 * R * n1);
        let escaped = 2.0 * R * incident * (direct * 400.0 + (tau - direct) * tw);
        let qwall = -2.0 * R * incident * (1.0 - direct) * (tw - 400.0);
        // Independent inversion of V*[a*(T1-T0)+b/2*(T1²-T0²)].
        let c0 = 1e6 + 1000.0 * tw;
        let tw1 = tw + (2.0 * qwall / 1e-6) / ((c0 * c0 + 2.0 * 1000.0 * qwall / 1e-6).sqrt() + c0);
        near(result.parcels()[&0].mass_kg(), n1 * M, 1e-18);
        near(result.parcels()[&0].state().state().temperature, tg, 3e-8);
        near(result.step().gas.audit.escaped_energy_j, escaped, 1e-11);
        near(result.step().wall_energy_change_j, qwall, 1e-9);
        near(result.walls()[&9].temperature_k(), tw1, 1e-10);
        assert!(result.step().energy_residual_j.abs() <= controls().exchange.energy_tolerance_j);
        let entropy = 1e-6 * (1e6 * ((tw1 - tw) / tw).ln_1p() + 1000.0 * (tw1 - tw));
        near(result.step().wall_entropy_change_j_k, entropy, 2e-12);
        // A hot wall loses entropy legitimately; vacuum has no invented entropy.
        assert_eq!(result.step().wall_entropy_change_j_k > 0.0, tw < 400.0);
        assert_eq!(result.walls()[&9].law().identity(), law.identity());
        assert_eq!(result.walls()[&9].mass_kg(), walls[&9].mass_kg());
        near(walls[&9].temperature_k(), tw, 1e-10);
    }
}

#[test]
fn accepted_finite_wall_temperature_drives_the_next_molecular_stage() {
    let source = material(1e6, 0.0, [200.0, 1800.0]);
    let law = FixedVolumeThermalLaw::try_new(&source, 1000.0, 200.0)
        .unwrap()
        .with_sensible_entropy();
    let walls = BTreeMap::from([(
        9,
        FixedVolumeChannelWall::try_new_isochoric(&law, 1e-8, 900.0).unwrap(),
    )]);
    let path = channel(2e-8);
    let gases = BTreeMap::from([(0, gas(400.0, 100_000.0, 0.0))]);
    let first = step_candidate(
        &gases,
        &walls,
        GasTransportPaths::default(),
        &[link(&path)],
        DT,
        &endpoints(&gases),
        controls(),
        || false,
    )
    .unwrap();
    let second = step_candidate(
        first.parcels(),
        first.walls(),
        GasTransportPaths::default(),
        &[link(&path)],
        DT,
        &endpoints(first.parcels()),
        controls(),
        || false,
    )
    .unwrap();
    let held = native::step_candidate(
        first.parcels(),
        GasTransportPaths::default(),
        &[VacuumChannelNetworkLink {
            id: 7,
            parcel: 0,
            channel: &path,
            area_m2: AREA,
            wall_temperature_k: 900.0,
        }],
        DT,
        &endpoints(first.parcels()),
        controls(),
        || false,
    )
    .unwrap();
    assert!(first.walls()[&9].temperature_k() < 900.0);
    assert!(second.walls()[&9].temperature_k() < first.walls()[&9].temperature_k());
    assert!(
        second.parcels()[&0].state().state().temperature
            < held.parcels()[&0].state().state().temperature
    );
    assert_eq!(
        second.step().gas.outlets[0].wall_temperature_k,
        first.walls()[&9].temperature_k()
    );
}

#[test]
fn shared_wall_pools_opposing_outlet_heat_before_inversion_in_canonical_order() {
    let source = material(1000.0, 0.0, [499.0, 501.0]);
    let law = FixedVolumeThermalLaw::try_new(&source, 1000.0, 500.0)
        .unwrap()
        .with_sensible_entropy();
    let walls = BTreeMap::from([(
        9,
        FixedVolumeChannelWall::try_new_isochoric(&law, 1e-6, 500.0).unwrap(),
    )]);
    let path = channel(1e-8);
    let gases = BTreeMap::from([
        (0, gas(800.0, 200_000.0, 0.0)),
        (1, gas(200.0, 100_000.0, 0.0)),
    ]);
    let a = link(&path);
    let b = FiniteVacuumChannelLink {
        id: 8,
        parcel: 1,
        ..a
    };
    let run = |links: &[FiniteVacuumChannelLink<'_>]| {
        step_candidate(
            &gases,
            &walls,
            GasTransportPaths::default(),
            links,
            DT,
            &endpoints(&gases),
            controls(),
            || false,
        )
        .unwrap()
    };
    let forward = run(&[a, b]);
    let reverse = run(&[b, a]);
    assert_eq!(forward.parcels(), reverse.parcels());
    assert_eq!(forward.step(), reverse.step());
    assert_eq!(forward.step().walls.len(), 1);
    near(forward.walls()[&9].temperature_k(), 500.0, 1e-8);
    for outlet in &forward.step().gas.outlets {
        assert!(
            law.temperature_k(outlet.wall.heat_into_wall_j / walls[&9].mass_kg())
                .is_err()
        );
    }
}

#[test]
fn source_reference_shifts_leave_finite_gas_and_wall_trajectories_unchanged() {
    let source = material(1e6, 0.0, [200.0, 1800.0]);
    let law = FixedVolumeThermalLaw::try_new(&source, 1000.0, 200.0)
        .unwrap()
        .with_sensible_entropy();
    let walls = BTreeMap::from([(
        9,
        FixedVolumeChannelWall::try_new_isochoric(&law, 1e-6, 900.0).unwrap(),
    )]);
    let mut previous = None;
    for datum in [0.0, 1e20] {
        let path = DiffuseMolecularChannel::try_new(
            interface(datum),
            CylindricalMolecularGeometry::try_new(1e-8, 2e-8, 1e-4, "synthetic channel").unwrap(),
        )
        .unwrap();
        let gases = BTreeMap::from([(0, gas(400.0, 100_000.0, datum))]);
        let result = step_candidate(
            &gases,
            &walls,
            GasTransportPaths::default(),
            &[link(&path)],
            DT,
            &endpoints(&gases),
            controls(),
            || false,
        )
        .unwrap();
        let values = [
            result.parcels()[&0].state().state().temperature,
            result.walls()[&9].temperature_k(),
            result.step().wall_energy_change_j,
        ];
        if let Some(old) = previous {
            for (value, original) in values.into_iter().zip(old) {
                near(value, original, 1e-10);
            }
        }
        previous = Some(values);
        assert!(result.step().energy_residual_j.abs() <= controls().exchange.energy_tolerance_j);
    }
}

#[test]
fn inactive_and_ballistic_outlets_do_not_query_an_unexposed_wall_domain() {
    let source = material(1e6, 0.0, [1800.0, 2000.0]);
    let law = FixedVolumeThermalLaw::try_new(&source, 1000.0, 1800.0)
        .unwrap()
        .with_sensible_entropy();
    let walls = BTreeMap::from([(
        9,
        FixedVolumeChannelWall::try_new_isochoric(&law, 1e-6, 1900.0).unwrap(),
    )]);
    let gases = BTreeMap::from([(0, gas(600.0, 150_000.0, 0.0))]);
    for (length, area) in [(0.0, AREA), (1e-8, 0.0)] {
        let path = channel(length);
        let outlet = FiniteVacuumChannelLink {
            area_m2: area,
            ..link(&path)
        };
        let result = step_candidate(
            &gases,
            &walls,
            GasTransportPaths::default(),
            &[outlet],
            DT,
            &endpoints(&gases),
            controls(),
            || false,
        )
        .unwrap();
        assert_eq!(
            result.walls()[&9].specific_energy_j_kg(),
            walls[&9].specific_energy_j_kg()
        );
        assert_eq!(result.step().wall_energy_change_j, 0.0);
        assert_eq!(result.step().wall_entropy_change_j_k, 0.0);
        if area == 0.0 {
            assert_eq!(result.parcels(), &gases);
        } else {
            assert!(result.parcels()[&0].mass_kg() < gases[&0].mass_kg());
        }
    }
}

#[test]
fn all_paths_share_gross_inventory_and_failures_leave_gas_and_wall_unchanged() {
    let source = material(1e6, 0.0, [200.0, 1800.0]);
    let law = FixedVolumeThermalLaw::try_new(&source, 1000.0, 200.0)
        .unwrap()
        .with_sensible_entropy();
    let walls = BTreeMap::from([(
        9,
        FixedVolumeChannelWall::try_new_isochoric(&law, 1e-6, 900.0).unwrap(),
    )]);
    let path = channel(2e-8);
    let gases = BTreeMap::from([
        (0, gas(400.0, 100_000.0, 0.0)),
        (1, gas(400.0, 100_000.0, 0.0)),
    ]);
    let original = gases.clone();
    let old_wall = walls[&9].specific_energy_j_kg();
    let native_paths = [EffusionNetworkLink {
        id: 99,
        parcels: [0, 1],
        interface: path.mouth(),
        area_m2: AREA,
    }];
    let fraction = incident_moles(400.0, 100_000.0) * M / gases[&0].mass_kg();
    let mut limited = controls();
    limited.exchange.maximum_species_fraction = 1.5 * fraction;
    assert!(matches!(
        step_candidate(
            &gases,
            &walls,
            GasTransportPaths {
                effusion: &native_paths,
                ..GasTransportPaths::default()
            },
            &[link(&path)],
            DT,
            &endpoints(&gases),
            limited,
            || false,
        ),
        Err(FiniteChannelError::Transport(GasTransportError::Exchange(
            SpeciesNetworkError::Withdrawal { .. }
        )))
    ));
    let tiny_wall = BTreeMap::from([(
        9,
        FixedVolumeChannelWall::try_new_isochoric(&law, 1e-11, 900.0).unwrap(),
    )]);
    assert!(matches!(
        step_candidate(
            &gases,
            &tiny_wall,
            GasTransportPaths::default(),
            &[link(&path)],
            DT,
            &endpoints(&gases),
            controls(),
            || false,
        ),
        Err(FiniteChannelError::WallDomain { .. } | FiniteChannelError::WallOvershoot { .. })
    ));
    let mut calls = 0;
    step_candidate(
        &gases,
        &walls,
        GasTransportPaths::default(),
        &[link(&path)],
        DT,
        &endpoints(&gases),
        controls(),
        || {
            calls += 1;
            false
        },
    )
    .unwrap();
    let mut cancelled_calls = 0;
    assert!(matches!(
        step_candidate(
            &gases,
            &walls,
            GasTransportPaths::default(),
            &[link(&path)],
            DT,
            &endpoints(&gases),
            controls(),
            || {
                cancelled_calls += 1;
                cancelled_calls == calls
            },
        ),
        Err(FiniteChannelError::Cancelled)
    ));
    assert_eq!(gases, original);
    assert_eq!(walls[&9].specific_energy_j_kg(), old_wall);
}
