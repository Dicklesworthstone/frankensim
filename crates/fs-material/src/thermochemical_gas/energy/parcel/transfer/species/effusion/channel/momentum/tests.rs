//! G1/G3: angular collimation, original molecular speeds, finite impulses and
//! complete stationary gas/channel momentum balances.
use super::*;
use super::super::{CylindricalMolecularGeometry, SpeciesThermalAccommodation};
use super::super::super::tests::{interface, pure};
use crate::thermochemical_gas::{
    ThermochemicalGasModel,
    energy::{GasEnergySolve, parcel::ClosedGasParcel},
    tests::{point, species_in_domain, transport},
};
use crate::thermochemical_gas::energy::parcel::transfer::species::network::{
    SpeciesNetworkSolve,
    transport::{
        GasTransportPaths, GasTransportSolve,
        effusion::channel::{
            self as finite,
            DiffuseChannelNetworkLink,
            vacuum::{self as vacuum, VacuumChannelNetworkLink},
        },
    },
};
use fs_thermochem::{Composition, CompositionBasis, mixture::FrozenIdealGasMixtureModelV1};
use std::collections::BTreeMap;

const AREA: f64 = 1e-10;
const D: f64 = 1e-9;
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}
fn base(gas: &ClosedGasParcel, h: f64) -> DiffuseMolecularChannel {
    DiffuseMolecularChannel::try_new(
        interface(gas.model().clone(), 1.0),
        CylindricalMolecularGeometry::try_new(
            D, D * h, 1e-7, "synthetic circular diffuse passage",
        ).unwrap(),
    ).unwrap()
}
fn momentum(channel: DiffuseMolecularChannel) -> DiffuseChannelMomentumInterface {
    DiffuseChannelMomentumInterface::try_new(
        channel, ClausingControls::default(),
        "synthetic diffuse angle and full-Maxwell/elastic endpoint speed laws",
        || false,
    ).unwrap()
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
fn closed(exterior: [f64; 2], wall: f64) {
    near(exterior[0] + exterior[1] + wall, 0.0, 1e-18);
}

#[test]
fn aperture_and_equilibrium_pressure_limits_preserve_original_molecular_moments() {
    let gas = pure(0.028, 2.5, 0.0, 0.0, 400.0, 100_000.0);
    let aperture = momentum(base(&gas, 0.0));
    let vacuum = aperture.sample_to_vacuum(&gas, AREA, f64::NAN, || false).unwrap();
    let thin = aperture.mouth.vacuum_flux(&vacuum.flux.transmitted, || false).unwrap();
    assert_eq!(vacuum.species[0].from_left.transmitted_n, thin.normal_momentum_rate_n());
    assert_eq!(vacuum.wall_force_n, 0.0);
    assert_eq!(vacuum.species[0].from_left.returned_n, 0.0);
    near(vacuum.exterior_forces_n[1], AREA * 100_000.0 / 2.0, 1e-19);
    closed(vacuum.exterior_forces_n, vacuum.wall_force_n);

    for h in [0.0, 1.0, 3.0] {
        let law = momentum(base(&gas, h));
        let pair = law.sample([&gas, &gas], AREA, 400.0, || false).unwrap();
        // Both molecular stresses remain when opposite mass fluxes cancel.
        assert_eq!(pair.flux.transmitted.net_mass_rate_left_to_right_kg_s, 0.0);
        near(pair.exterior_forces_n[0], -AREA * 100_000.0, 2e-11);
        near(pair.exterior_forces_n[1], AREA * 100_000.0, 2e-11);
        assert_eq!(pair.wall_force_n, 0.0);
        assert!(pair.species[0].from_left.transmitted_n > 0.0);
        assert_eq!(pair.species[0].from_left, pair.species[0].from_right);
        closed(pair.exterior_forces_n, pair.wall_force_n);

        let zero = law.sample_to_vacuum(&gas, 0.0, f64::NAN, || false).unwrap();
        assert!(zero.species.is_empty());
        assert_eq!(zero.exterior_forces_n, [0.0; 2]);
        assert_eq!(zero.wall_force_n, 0.0);
    }
}

#[test]
fn channel_collimation_and_wall_speed_change_thrust_without_changing_mass() {
    let gas = pure(0.028, 3.5, 0.001, 0.0, 400.0, 100_000.0);
    let law = momentum(base(&gas, 2.0));
    let at_gas = law.sample_to_vacuum(&gas, AREA, 400.0, || false).unwrap();
    let hot = law.sample_to_vacuum(&gas, AREA, 900.0, || false).unwrap();
    let cold = law.sample_to_vacuum(&gas, AREA, 300.0, || false).unwrap();
    let tau = law.channel().geometry().transmission_probability();
    let a = at_gas.species[0].from_left;
    let b = hot.species[0].from_left;
    let c = cold.species[0].from_left;
    assert_eq!(at_gas.flux.transmitted, hot.flux.transmitted);
    assert_eq!(at_gas.flux.transmitted, cold.flux.transmitted);
    near(a.incident_n, AREA * 100_000.0 / 2.0, 1e-19);
    assert!(a.transmitted_n > tau * a.incident_n,
        "a finite bore collimates; conductance alone is not thrust");
    near(a.transmitted_n + a.returned_n, a.incident_n, 2e-11);
    near(b.returned_n / a.returned_n, (900.0_f64 / 400.0).sqrt(), 1e-13);
    near(c.returned_n / a.returned_n, (300.0_f64 / 400.0).sqrt(), 1e-13);
    // The direct contribution never encounters the wall. Only the residual
    // outgoing contribution follows the Maxwellian sqrt(T_wall) speed scale.
    let direct = a.incident_n * 1.5 * law.moments().direct_axial_cosine();
    near((b.transmitted_n - direct) / (a.transmitted_n - direct), 1.5, 1e-12);
    assert!(hot.flux.wall.heat_into_escaped_w > 0.0);
    for result in [&at_gas, &hot, &cold] {
        closed(result.exterior_forces_n, result.wall_force_n);
        near(result.exterior_forces_n[0] + result.wall_force_n,
            -result.exterior_forces_n[1], 1e-18);
    }
    let reversed = law.sample([&gas, &gas], AREA, 900.0, || false).unwrap();
    assert_eq!(reversed.wall_force_n, 0.0);
}

#[test]
fn species_speed_laws_keep_elastic_sources_inactive_and_refuse_unresolved_histories() {
    let model = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![
                species_in_domain("A", 0.028, 2.5, 0.0, 0.0, [200.0, 500.0]),
                species_in_domain("B", 0.004, 2.5, 0.0, 1e24, [200.0, 2000.0]),
            ],
            Composition::new(CompositionBasis::MoleFraction, vec![0.25, 0.75]).unwrap(),
        ).unwrap(),
        vec![transport("A", 1.7e-5, 0.025), transport("B", 2e-5, 0.15)],
        None,
    ).unwrap();
    let gas = ClosedGasParcel::new(model, point(400.0, 100_000.0), 0.001).unwrap();
    let table = base(&gas, 1.0).with_species_thermal_accommodation(
        vec![
            SpeciesThermalAccommodation::try_new(
                SpeciesId::new("B").unwrap(), 1.0, "full Maxwellian species B",
            ).unwrap(),
            SpeciesThermalAccommodation::try_new(
                SpeciesId::new("A").unwrap(), 0.0, "speed-preserving elastic species A",
            ).unwrap(),
        ],
        ClausingControls::default(), || false,
    ).unwrap();
    let law = momentum(table);
    let a = law.sample_to_vacuum(&gas, AREA, 400.0, || false).unwrap();
    let b = law.sample_to_vacuum(&gas, AREA, 900.0, || false).unwrap();
    assert_eq!(a.species[0].from_left, b.species[0].from_left);
    near(b.species[1].from_left.returned_n / a.species[1].from_left.returned_n, 1.5, 1e-13);
    assert_eq!(a.flux.transmitted, b.flux.transmitted);
    let geometry = law.channel().geometry();
    for row in law.channel().species_thermal_accommodation().unwrap() {
        if row.qualification().coefficient() == 1.0 {
            assert_eq!(row.accommodation().reflected_heat_probability(), geometry.reflected_probability());
            assert_eq!(row.accommodation().transmitted_heat_probability(), geometry.wall_transmitted_probability());
        } else {
            assert_eq!(row.accommodation().reflected_heat_probability(), 0.0);
        }
    }

    let pure_gas = pure(0.028, 2.5, 0.0, 0.0, 400.0, 100_000.0);
    let elastic = momentum(base(&pure_gas, 1.0).with_thermal_accommodation(
        0.0, ClausingControls::default(), "elastic diffuse speed law", || false,
    ).unwrap());
    let inactive = elastic.sample_to_vacuum(&pure_gas, AREA, f64::NAN, || false).unwrap();
    assert_eq!(inactive.flux.wall.heat_into_wall_w, 0.0);
    assert!(inactive.species[0].from_left.returned_n > 0.0);
    for h in [0.0, 1.0] {
        let partial = base(&pure_gas, h).with_thermal_accommodation(
            0.5, ClausingControls::default(), "Bernoulli collision energy law", || false,
        ).unwrap();
        let result = DiffuseChannelMomentumInterface::try_new(
            partial, ClausingControls::default(), "explicit speed qualification", || false,
        );
        assert_eq!(result.is_ok(), h == 0.0);
    }
}

#[test]
fn accepted_vacuum_and_reciprocal_receipts_reuse_actual_species_amounts() {
    let gas = pure(0.028, 2.5, 0.0, 0.0, 400.0, 100_000.0);
    let law = momentum(base(&gas, 1.0));
    let gases = BTreeMap::from([(3, gas.clone())]);
    let endpoint = BTreeMap::from([(3, energy())]);
    let duration = 0.7;
    let accepted = vacuum::step_candidate(
        &gases, GasTransportPaths::default(),
        &[VacuumChannelNetworkLink {
            id: 9, parcel: 3, channel: law.channel(), area_m2: AREA, wall_temperature_k: 900.0,
        }],
        duration, &endpoint, controls(), || false,
    ).unwrap();
    let outlet = &accepted.step().outlets[0];
    let p = law.vacuum_outlet(outlet, || false).unwrap();
    let rates = law.sample_to_vacuum(&gas, AREA, 900.0, || false).unwrap();
    near(p.exterior_impulses_n_s[1], duration * rates.exterior_forces_n[1], 1e-18);
    near(p.wall_impulse_n_s, duration * rates.wall_force_n, 1e-18);
    closed(p.exterior_impulses_n_s, p.wall_impulse_n_s);
    let incoming_speed = (core::f64::consts::PI * R * 400.0 / (2.0 * 0.028)).sqrt();
    near(p.species[0].from_left.incident_n_s,
        (outlet.escaped_species[0].mass_kg + outlet.wall.reflected_species[0].mass_kg)
            * incoming_speed, 1e-18);
    let original = outlet.clone();
    let mut duplicate = outlet.clone();
    duplicate.escaped_species.push(duplicate.escaped_species[0].clone());
    assert!(law.vacuum_outlet(&duplicate, || false).is_err());
    let mut missing = outlet.clone();
    missing.escaped_species.clear();
    assert!(law.vacuum_outlet(&missing, || false).is_err());
    let mut foreign = outlet.clone();
    foreign.wall.reflected_species[0].species = SpeciesId::new("foreign").unwrap();
    assert!(law.vacuum_outlet(&foreign, || false).is_err());
    let mut unqualified = outlet.clone();
    unqualified.channel = base(&gas, 1.0);
    assert!(law.vacuum_outlet(&unqualified, || false).is_err());
    assert_eq!(outlet, &original);

    let right = pure(0.028, 2.5, 0.0, 0.0, 600.0, 70_000.0);
    let pair = BTreeMap::from([(3, gas.clone()), (4, right.clone())]);
    let pair_step = finite::step_candidate(
        &pair, GasTransportPaths::default(),
        &[DiffuseChannelNetworkLink {
            id: 8, parcels: [3, 4], channel: law.channel(), area_m2: AREA,
            wall_temperature_k: 900.0,
        }],
        duration, &BTreeMap::from([(3, energy()), (4, energy())]), controls(), || false,
    ).unwrap();
    let finite_p = law.reciprocal_step(pair_step.step(), 8, || false).unwrap();
    let pair_rates = law.sample([&gas, &right], AREA, 900.0, || false).unwrap();
    for side in 0..2 {
        near(finite_p.exterior_impulses_n_s[side], duration * pair_rates.exterior_forces_n[side], 1e-18);
    }
    closed(finite_p.exterior_impulses_n_s, finite_p.wall_impulse_n_s);

    let mut calls = 0;
    law.vacuum_outlet(outlet, || { calls += 1; false }).unwrap();
    let mut seen = 0;
    assert!(matches!(law.vacuum_outlet(outlet, || {
        seen += 1; seen == calls
    }), Err(GasEnergyError::Cancelled)));
    assert_eq!(outlet, &original);
}
