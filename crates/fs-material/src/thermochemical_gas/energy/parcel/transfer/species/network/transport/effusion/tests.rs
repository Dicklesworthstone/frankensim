//! G1/G3 ballistic finite-reservoir oracles through the original shared owners.
use super::*;
use crate::thermochemical_gas::{
    ThermochemicalGasModel,
    energy::{
        GasEnergySolve,
        parcel::transfer::{
            nozzle::{GasNozzle, network::snapshot::SnapshotNozzleLink},
            species::{
                effusion::{EffusiveRegime, EffusiveSpeciesTransmission, EffusiveSupport},
                network::{
                    SpeciesNetworkError, SpeciesNetworkSolve,
                    transport::{
                        GasTransportError, GasTransportPaths, GasTransportSolve,
                        step_candidate_with_paths,
                    },
                },
            },
        },
    },
    tests::{point, species_with_energy_reference, transport_in_domain},
};
use fs_thermochem::{
    Composition, CompositionBasis, SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
    mixture::FrozenIdealGasMixtureModelV1,
};

const M: f64 = 0.028;
const VOLUME: f64 = 0.001;

fn model(datum: f64) -> ThermochemicalGasModel {
    ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![species_with_energy_reference("A", M, 2.5, 0.0, datum)],
            Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
        )
        .unwrap(),
        vec![transport_in_domain(
            "A",
            1.7e-5,
            0.025,
            [200.0, 2000.0],
            [1.0, 1e7],
        )],
        None,
    )
    .unwrap()
}
fn gas(t: f64, p: f64, datum: f64) -> ClosedGasParcel {
    ClosedGasParcel::new(model(datum), point(t, p), VOLUME).unwrap()
}
fn interface(datum: f64) -> EffusiveInterface {
    EffusiveInterface::try_new(
        model(datum),
        vec![
            EffusiveSpeciesTransmission::try_new(
                SpeciesId::new("A").unwrap(),
                1.0,
                "synthetic reciprocal aperture",
            )
            .unwrap(),
        ],
        EffusiveSupport {
            temperature_k: (200.0, 1800.0),
            pressure_pa: (1.0, 1e7),
        },
        EffusiveRegime::try_new(
            1e-8,
            1e-11,
            1e-4,
            "synthetic collisionless thin-array assumption; no measured regime authority",
        )
        .unwrap(),
    )
    .unwrap()
}
fn energy() -> GasEnergySolve {
    GasEnergySolve {
        minimum_temperature_k: 200.0,
        maximum_temperature_k: 1800.0,
        energy_tolerance_j_kg: 1e-7,
        maximum_iterations: 128,
    }
}
fn controls() -> GasTransportSolve {
    GasTransportSolve {
        exchange: SpeciesNetworkSolve {
            maximum_species_fraction: 0.1,
            mass_tolerance_kg: 1e-13,
            energy_tolerance_j: 1e-8,
        },
        entropy_tolerance_j_k: 1e-10,
    }
}
fn link(interface: &EffusiveInterface, area_m2: f64) -> EffusionNetworkLink<'_> {
    EffusionNetworkLink {
        id: 8,
        parcels: [0, 1],
        interface,
        area_m2,
    }
}
fn rate(t: f64, p: f64, area: f64) -> f64 {
    area * p / (2.0 * core::f64::consts::PI * M * R * t).sqrt()
}
fn near(actual: f64, expected: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute,
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn monatomic_finite_update_matches_two_independent_ballistic_energy_beams() {
    let parcels = BTreeMap::from([(0, gas(600.0, 150.0, 0.0)), (1, gas(300.0, 100.0, 0.0))]);
    let original = parcels.clone();
    let surface = interface(0.0);
    let (area, dt) = (1e-8, 1.0);
    let paths = [link(&surface, area)];
    let endpoints = BTreeMap::from([(0, energy()), (1, energy())]);
    let candidate = step_candidate_with_paths(
        &parcels,
        GasTransportPaths {
            effusion: &paths,
            ..Default::default()
        },
        dt,
        &endpoints,
        controls(),
        || false,
    )
    .unwrap();
    let (nl, nr) = (parcels[&0].mass_kg() / M, parcels[&1].mass_kg() / M);
    let (forward, reverse) = (dt * rate(600.0, 150.0, area), dt * rate(300.0, 100.0, area));
    // Independent monatomic result: each crossing particle carries 2RT;
    // stored internal energy uses Cv=3R/2, not enthalpy or a second heat source.
    let energy_lr = 2.0 * R * (forward * 600.0 - reverse * 300.0);
    let left_t = (1.5 * R * nl * 600.0 - energy_lr) / (1.5 * R * (nl - forward + reverse));
    let right_t = (1.5 * R * nr * 300.0 + energy_lr) / (1.5 * R * (nr + forward - reverse));
    near(
        candidate.parcels()[&0].mass_kg(),
        M * (nl - forward + reverse),
        1e-18,
    );
    near(
        candidate.parcels()[&1].mass_kg(),
        M * (nr + forward - reverse),
        1e-18,
    );
    near(
        candidate.parcels()[&0].state().state().temperature,
        left_t,
        2e-8,
    );
    near(
        candidate.parcels()[&1].state().state().temperature,
        right_t,
        2e-8,
    );
    near(
        candidate.step().exchange.parcels[0].energy_change_j,
        -energy_lr,
        2e-10,
    );
    near(
        candidate.step().exchange.parcels[1].energy_change_j,
        energy_lr,
        2e-10,
    );
    assert_eq!(
        candidate.step().exchange.fluxes.len(),
        2,
        "both original beams survive"
    );
    assert!(candidate.step().entropy_change_j_k > 0.0);
    assert!(
        candidate.step().exchange.total_energy_residual_j.abs()
            <= controls().exchange.energy_tolerance_j
    );
    assert_eq!(parcels, original);
}

#[test]
fn zero_net_mass_still_transfers_heat_and_survives_large_formation_reference() {
    let mut temperatures = Vec::new();
    for datum in [0.0, 1e12] {
        // sqrt(1200/300)=2 and pressure ratio2 give exact transpiration balance.
        let parcels = BTreeMap::from([
            (0, gas(1200.0, 200.0, datum)),
            (1, gas(300.0, 100.0, datum)),
        ]);
        let surface = interface(datum);
        let paths = [link(&surface, 1e-8)];
        let endpoints = BTreeMap::from([(0, energy()), (1, energy())]);
        let result = step_candidate_with_paths(
            &parcels,
            GasTransportPaths {
                effusion: &paths,
                ..Default::default()
            },
            1.0,
            &endpoints,
            controls(),
            || false,
        )
        .unwrap();
        let heat = rate(300.0, 100.0, 1e-8) * 2.0 * R * (1200.0 - 300.0);
        near(result.parcels()[&0].mass_kg(), parcels[&0].mass_kg(), 1e-19);
        near(result.parcels()[&1].mass_kg(), parcels[&1].mass_kg(), 1e-19);
        near(
            result.step().exchange.parcels[0].energy_change_j,
            -heat,
            2e-9,
        );
        near(
            result.step().exchange.parcels[1].energy_change_j,
            heat,
            2e-9,
        );
        assert!(result.step().exchange.parcels[0].net_additional_heat_j > 0.0);
        assert!(result.step().entropy_change_j_k > 0.0);
        let pair = [
            result.parcels()[&0].state().state().temperature,
            result.parcels()[&1].state().state().temperature,
        ];
        assert!(pair[0] < 1200.0 && pair[1] > 300.0);
        temperatures.push(pair);
    }
    for side in 0..2 {
        near(temperatures[0][side], temperatures[1][side], 2e-8);
    }
}

#[test]
fn thermal_transpiration_can_move_mass_toward_higher_total_pressure() {
    let parcels = BTreeMap::from([(0, gas(1200.0, 150.0, 0.0)), (1, gas(300.0, 100.0, 0.0))]);
    let surface = interface(0.0);
    let paths = [link(&surface, 1e-8)];
    let endpoints = BTreeMap::from([(0, energy()), (1, energy())]);
    let result = step_candidate_with_paths(
        &parcels,
        GasTransportPaths {
            effusion: &paths,
            ..Default::default()
        },
        0.1,
        &endpoints,
        controls(),
        || false,
    )
    .unwrap();
    assert!(
        result.step().effusion[0]
            .flux
            .net_mass_rate_left_to_right_kg_s
            < 0.0
    );
    assert!(result.parcels()[&0].mass_kg() > parcels[&0].mass_kg());
    assert!(result.parcels()[&1].mass_kg() < parcels[&1].mass_kg());
    result
        .step()
        .admit_endpoint_driving_forces(result.parcels(), || false)
        .unwrap();
    assert!(result.step().entropy_change_j_k > 0.0);
}

#[test]
fn effusive_and_nozzle_beams_share_the_same_original_gross_allowance() {
    let parcels = BTreeMap::from([
        (0, gas(400.0, 200_000.0, 0.0)),
        (1, gas(400.0, 100_000.0, 0.0)),
        (2, gas(400.0, 100_000.0, 0.0)),
    ]);
    let endpoints: BTreeMap<_, _> = parcels.keys().map(|&id| (id, energy())).collect();
    let surface = interface(0.0);
    let dt = 0.1;
    let target = 0.06 * parcels[&0].mass_kg();
    let area = target / (dt * M * rate(400.0, 200_000.0, 1.0));
    let effusive = [link(&surface, area)];
    let mut nozzle = SnapshotNozzleLink {
        id: 5,
        parcels: [0, 2],
        nozzle: GasNozzle {
            area_m2: 1e-8,
            discharge_coefficient: 0.8,
        },
        nozzle_solve: energy(),
    };
    let sampled = parcels[&0]
        .nozzle_flow_to(100_000.0, nozzle.nozzle, energy(), || false)
        .unwrap();
    nozzle.nozzle.area_m2 *= target / (sampled.mass_flow_kg_s * dt);
    step_candidate_with_paths(
        &parcels,
        GasTransportPaths {
            effusion: &effusive,
            ..Default::default()
        },
        dt,
        &endpoints,
        controls(),
        || false,
    )
    .unwrap();
    step_candidate_with_paths(
        &parcels,
        GasTransportPaths {
            nozzles: &[nozzle],
            ..Default::default()
        },
        dt,
        &endpoints,
        controls(),
        || false,
    )
    .unwrap();
    assert!(matches!(
        step_candidate_with_paths(
            &parcels,
            GasTransportPaths {
                nozzles: &[nozzle],
                effusion: &effusive,
                ..Default::default()
            },
            dt,
            &endpoints,
            controls(),
            || false
        ),
        Err(GasTransportError::Exchange(
            SpeciesNetworkError::Withdrawal { id: 0, .. }
        ))
    ));
}

#[test]
fn equilibrium_closed_paths_ordering_and_late_cancellation_are_atomic() {
    let parcels = BTreeMap::from([(0, gas(400.0, 100.0, 0.0)), (1, gas(400.0, 100.0, 0.0))]);
    let original = parcels.clone();
    let surface = interface(0.0);
    let endpoints = BTreeMap::from([(0, energy()), (1, energy())]);
    let paths = [link(&surface, 1e-8)];
    let equilibrium = step_candidate_with_paths(
        &parcels,
        GasTransportPaths {
            effusion: &paths,
            ..Default::default()
        },
        1.0,
        &endpoints,
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(equilibrium.parcels(), &parcels);
    assert_eq!(equilibrium.step().exchange.fluxes.len(), 2);
    assert_eq!(equilibrium.step().entropy_change_j_k, 0.0);
    let closed = [link(&surface, 0.0)];
    let result = step_candidate_with_paths(
        &parcels,
        GasTransportPaths {
            effusion: &closed,
            ..Default::default()
        },
        1.0,
        &endpoints,
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(result.parcels(), &parcels);
    assert!(result.step().exchange.fluxes.is_empty());
    assert!(result.step().exchange.heat.is_empty());
    let mut a = link(&surface, 3e-9);
    a.id = 1;
    let b = link(&surface, 7e-9);
    let run = |links: &[EffusionNetworkLink<'_>]| {
        step_candidate_with_paths(
            &parcels,
            GasTransportPaths {
                effusion: links,
                ..Default::default()
            },
            1.0,
            &endpoints,
            controls(),
            || false,
        )
        .unwrap()
    };
    assert_eq!(run(&[a, b]), run(&[b, a]));
    let mut polls = 0;
    step_candidate_with_paths(
        &parcels,
        GasTransportPaths {
            effusion: &paths,
            ..Default::default()
        },
        1.0,
        &endpoints,
        controls(),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let mut count = 0;
    assert_eq!(
        step_candidate_with_paths(
            &parcels,
            GasTransportPaths {
                effusion: &paths,
                ..Default::default()
            },
            1.0,
            &endpoints,
            controls(),
            || {
                count += 1;
                count == polls
            }
        ),
        Err(GasTransportError::Cancelled)
    );
    assert_eq!(parcels, original);
}
