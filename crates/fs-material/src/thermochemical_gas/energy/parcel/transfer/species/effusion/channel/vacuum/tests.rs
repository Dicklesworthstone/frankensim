//! Independent molecule-by-molecule vacuum-channel balances.

use super::super::super::{
    R,
    collision::{EffusiveCollisionModel, EffusiveCollisionPair},
    tests::{interface, near, pure},
};
use super::super::CylindricalMolecularGeometry;
use super::*;
use crate::thermochemical_gas::tests::{model, point};

const AREA: f64 = 1e-10;
const DIAMETER: f64 = 1e-9;

fn channel(gas: &ClosedGasParcel, length_m: f64) -> DiffuseMolecularChannel {
    DiffuseMolecularChannel::try_new(
        interface(gas.model().clone(), 1.0),
        CylindricalMolecularGeometry::try_new(
            DIAMETER,
            length_m,
            1e-7,
            "synthetic vacuum channel qualification",
        )
        .unwrap(),
    )
    .unwrap()
}

fn incident_moles(pressure: f64, temperature: f64, molar_mass: f64) -> f64 {
    AREA * pressure / (2.0 * core::f64::consts::PI * molar_mass * R * temperature).sqrt()
}

#[test]
fn returned_gas_escaped_molecules_and_wall_close_an_independent_energy_ledger() {
    let mass = 0.004;
    let temperature = 400.0;
    let pressure = 100_000.0;
    let gas = pure(mass, 2.5, 0.0, 0.0, temperature, pressure);
    let original = gas.clone();
    let law = channel(&gas, DIAMETER);
    let q = incident_moles(pressure, temperature, mass);
    let tau = 18.0 / 35.0;
    let direct = 0.171_572_875_253_809_9;
    let source_energy = 2.0 * R * temperature;
    for wall_temperature in [300.0, temperature, 900.0] {
        let flux = law
            .sample_to_vacuum(&gas, AREA, wall_temperature, || false)
            .unwrap();
        let wall_energy = 2.0 * R * wall_temperature;
        let change = wall_energy - source_energy;
        assert_eq!(flux.transmitted.species.len(), 1);
        assert_eq!(flux.wall.reflected_species.len(), 1);
        near(flux.transmitted.molar_rate_mol_s, tau * q);
        near(flux.transmitted.mass_rate_kg_s, tau * q * mass);
        near(
            flux.wall.reflected_species[0].mass_rate_kg_s,
            (1.0 - tau) * q * mass,
        );
        near(flux.wall.heat_into_gas_w, (1.0 - tau) * q * change);
        near(flux.wall.heat_into_escaped_w, (tau - direct) * q * change);
        near(
            flux.escaped_energy_rate_w,
            q * (direct * source_energy + (tau - direct) * wall_energy),
        );
        let gas_power = -flux.transmitted.energy_rate_w + flux.wall.heat_into_gas_w;
        near(gas_power, q * (-source_energy + (1.0 - tau) * wall_energy));
        near(flux.wall.heat_into_wall_w, -(1.0 - direct) * q * change);
        let residual = gas_power + flux.escaped_energy_rate_w + flux.wall.heat_into_wall_w;
        let scale =
            gas_power.abs() + flux.escaped_energy_rate_w.abs() + flux.wall.heat_into_wall_w.abs();
        assert!(residual.abs() <= 8e-12 * scale);
        if wall_temperature == temperature {
            assert_eq!(flux.wall.heat_into_gas_w, 0.0);
            assert_eq!(flux.wall.heat_into_escaped_w, 0.0);
            assert_eq!(flux.wall.heat_into_wall_w, 0.0);
        }
    }
    assert_eq!(gas, original);
}

#[test]
fn nonlinear_internal_modes_and_large_formation_references_preserve_wall_heat() {
    let mass = 0.028;
    let temperature = 400.0;
    let wall_temperature = 900.0;
    let pressure = 100_000.0;
    let cp_r = 3.5;
    let slope = 0.001;
    let q = incident_moles(pressure, temperature, mass);
    let change = R
        * (wall_temperature - temperature)
        * (cp_r - 0.5 + 0.5 * slope * (wall_temperature + temperature));
    let tau = 18.0 / 35.0;
    let direct = 0.171_572_875_253_809_9;
    let mut baseline = None;
    for reference in [0.0, 1e24, -1e24] {
        let gas = pure(mass, cp_r, slope, reference, temperature, pressure);
        let law = channel(&gas, DIAMETER);
        let flux = law
            .sample_to_vacuum(&gas, AREA, wall_temperature, || false)
            .unwrap();
        near(flux.wall.heat_into_gas_w, (1.0 - tau) * q * change);
        near(flux.wall.heat_into_escaped_w, (tau - direct) * q * change);
        near(flux.wall.heat_into_wall_w, -(1.0 - direct) * q * change);
        let heat = [
            flux.wall.heat_into_gas_w,
            flux.wall.heat_into_escaped_w,
            flux.wall.heat_into_wall_w,
        ];
        if let Some(previous) = baseline {
            for (actual, expected) in heat.into_iter().zip(previous) {
                near(actual, expected);
            }
        } else {
            baseline = Some(heat);
        }
        if reference < 0.0 {
            assert!(flux.escaped_energy_rate_w < 0.0);
        }
    }
}

#[test]
fn two_species_preserve_their_graham_rates_and_distinct_caloric_owners() {
    let gas =
        ClosedGasParcel::new(model(vec![0.75, 0.25], 0.0), point(400.0, 100_000.0), 0.001).unwrap();
    let law = channel(&gas, DIAMETER);
    let wall_temperature = 900.0;
    let flux = law
        .sample_to_vacuum(&gas, AREA, wall_temperature, || false)
        .unwrap();
    let tau = 18.0 / 35.0;
    let direct = 0.171_572_875_253_809_9;
    let mut expected_gas = 0.0;
    let mut expected_escaped = 0.0;
    let mut expected_energy = 0.0;
    for (index, (mass, fraction, cp_r)) in [(0.028, 0.75, 3.5), (0.004, 0.25, 2.5)]
        .into_iter()
        .enumerate()
    {
        let q = incident_moles(fraction * 100_000.0, 400.0, mass);
        let source_energy = R * (cp_r - 0.5) * 400.0;
        let wall_energy = R * (cp_r - 0.5) * wall_temperature;
        near(
            flux.transmitted.species[index].outward.molar_rate_mol_s,
            tau * q,
        );
        near(
            flux.wall.reflected_species[index].mass_rate_kg_s,
            (1.0 - tau) * q * mass,
        );
        expected_gas += (1.0 - tau) * q * (wall_energy - source_energy);
        expected_escaped += (tau - direct) * q * (wall_energy - source_energy);
        expected_energy += q * (direct * source_energy + (tau - direct) * wall_energy);
    }
    near(flux.wall.heat_into_gas_w, expected_gas);
    near(flux.wall.heat_into_escaped_w, expected_escaped);
    near(flux.escaped_energy_rate_w, expected_energy);
    near(flux.wall.heat_into_wall_w, -expected_gas - expected_escaped);
}

#[test]
fn one_way_accommodation_is_shared_with_the_two_reservoir_channel() {
    let gas = pure(0.028, 3.5, 0.001, 0.0, 400.0, 100_000.0);
    let law = channel(&gas, DIAMETER);
    let vacuum = law.sample_to_vacuum(&gas, AREA, 900.0, || false).unwrap();
    let pair = law.sample(&gas, &gas, AREA, 900.0, || false).unwrap();
    assert_eq!(
        vacuum.transmitted.species[0].outward,
        pair.transmitted.species[0].left_to_right
    );
    near(
        pair.wall.reflected_species[0].mass_rate_kg_s[0],
        vacuum.wall.reflected_species[0].mass_rate_kg_s,
    );
    near(
        pair.wall.heat_into_gas_w[0],
        vacuum.wall.heat_into_gas_w + vacuum.wall.heat_into_escaped_w,
    );
    near(
        pair.wall.heat_into_wall_w,
        2.0 * vacuum.wall.heat_into_wall_w,
    );
}

#[test]
fn zero_length_recovers_the_original_vacuum_aperture_without_a_wall_query() {
    let gas = pure(0.028, 3.5, 0.001, 0.0, 400.0, 100_000.0);
    let law = channel(&gas, 0.0);
    let original = gas
        .effusive_flux_to_vacuum(law.mouth(), AREA, || false)
        .unwrap();
    let flux = law
        .sample_to_vacuum(&gas, AREA, f64::NAN, || false)
        .unwrap();
    assert_eq!(flux.transmitted, original);
    assert_eq!(flux.escaped_energy_rate_w, original.energy_rate_w);
    assert!(flux.wall.reflected_species.is_empty());
    assert_eq!(flux.wall.heat_into_gas_w, 0.0);
    assert_eq!(flux.wall.heat_into_escaped_w, 0.0);
    assert_eq!(flux.wall.heat_into_wall_w, 0.0);
}

#[test]
fn closed_area_skips_foreign_gas_and_wall_but_active_outlets_admit_both() {
    let gas = pure(0.028, 3.5, 0.0, 0.0, 400.0, 100_000.0);
    let foreign = pure(0.004, 2.5, 0.0, 0.0, 1000.0, 100_000.0);
    let law = channel(&gas, DIAMETER);
    let closed = law
        .sample_to_vacuum(&foreign, 0.0, f64::NAN, || false)
        .unwrap();
    assert!(closed.transmitted.species.is_empty());
    assert!(closed.wall.reflected_species.is_empty());
    assert_eq!(closed.escaped_energy_rate_w, 0.0);
    assert_eq!(closed.wall.heat_into_gas_w, 0.0);
    assert_eq!(closed.wall.heat_into_escaped_w, 0.0);
    assert_eq!(closed.wall.heat_into_wall_w, 0.0);
    assert!(
        law.sample_to_vacuum(&foreign, AREA, 600.0, || false)
            .is_err()
    );
    for bad_area in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            law.sample_to_vacuum(&gas, bad_area, 600.0, || false)
                .is_err()
        );
    }
    for bad_temperature in [199.0, 2001.0, f64::NAN, f64::INFINITY] {
        assert!(
            law.sample_to_vacuum(&gas, AREA, bad_temperature, || false)
                .is_err()
        );
    }
    assert!(matches!(
        law.sample_to_vacuum(&gas, 0.0, f64::NAN, || true),
        Err(GasEnergyError::Cancelled)
    ));
    let mut polls = 0;
    law.sample_to_vacuum(&gas, AREA, 900.0, || {
        polls += 1;
        false
    })
    .unwrap();
    let mut seen = 0;
    assert!(matches!(
        law.sample_to_vacuum(&gas, AREA, 900.0, || {
            seen += 1;
            seen == polls
        }),
        Err(GasEnergyError::Cancelled)
    ));
}

#[test]
fn integrated_wall_energy_uses_the_same_rounded_mass_as_the_vacuum_boundary() {
    let gas = pure(0.028, 3.5, 0.001, 0.0, 400.0, 100_000.0);
    let law = channel(&gas, DIAMETER);
    let flux = law.sample_to_vacuum(&gas, AREA, 900.0, || false).unwrap();
    let duration = 0.123_456_789;
    let exchange = law
        .integrate_vacuum_wall_exchange(&gas, &flux.transmitted, 900.0, duration, || false)
        .unwrap();
    let escaped_mass = positive_product_ratio(
        &[flux.transmitted.species[0].outward.mass_rate_kg_s, duration],
        &[],
    )
    .unwrap();
    let tau = 18.0 / 35.0;
    let direct = 0.171_572_875_253_809_9;
    let change = (R / 0.028) * 500.0 * (3.0 + 0.0005 * 1300.0);
    near(
        exchange.reflected_species[0].mass_kg,
        escaped_mass * (1.0 - tau) / tau,
    );
    near(
        exchange.heat_into_gas_j,
        escaped_mass * (1.0 - tau) / tau * change,
    );
    near(
        exchange.heat_into_escaped_j,
        escaped_mass * (tau - direct) / tau * change,
    );
    near(
        exchange.heat_into_wall_j,
        -exchange.heat_into_gas_j - exchange.heat_into_escaped_j,
    );
    for invalid_duration in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            law.integrate_vacuum_wall_exchange(
                &gas,
                &flux.transmitted,
                900.0,
                invalid_duration,
                || false
            )
            .is_err()
        );
    }
}

#[test]
fn vacuum_does_not_disable_the_emitters_full_channel_collision_length() {
    let gas = pure(0.028, 3.5, 0.0, 0.0, 400.0, 100_000.0);
    let mouth = interface(gas.model().clone(), 1.0);
    let id = gas.model().thermodynamics().components()[0]
        .species()
        .clone();
    // At this source state lambda is about 49 nm: the 1 nm mouth admits,
    // but a 10 nm channel refuses its full collisional length.
    let collision = EffusiveCollisionModel::try_new(
        gas.model().clone(),
        vec![
            EffusiveCollisionPair::try_new(
                [id.clone(), id],
                8e-19,
                "synthetic dense collision partners",
            )
            .unwrap(),
        ],
        mouth.support(),
        "synthetic actual-state collision qualification",
    )
    .unwrap();
    let mouth = mouth.with_collision_model(collision).unwrap();
    mouth.admit_parcel(&gas, || false).unwrap();
    let law = DiffuseMolecularChannel::try_new(
        mouth,
        CylindricalMolecularGeometry::try_new(
            DIAMETER,
            10.0 * DIAMETER,
            1e-7,
            "synthetic geometry qualification",
        )
        .unwrap(),
    )
    .unwrap();
    assert!(law.sample_to_vacuum(&gas, AREA, 900.0, || false).is_err());
    assert!(law.sample_to_vacuum(&gas, 0.0, f64::NAN, || false).is_ok());
}

#[test]
fn a_source_subset_does_not_invent_absent_species_or_a_vacuum_state() {
    let mixed =
        ClosedGasParcel::new(model(vec![0.75, 0.25], 0.0), point(400.0, 100_000.0), 0.001).unwrap();
    let law = channel(&mixed, DIAMETER);
    let pure_b = mixed
        .model()
        .with_composition(
            fs_thermochem::Composition::new(
                fs_thermochem::CompositionBasis::MoleFraction,
                vec![0.0, 1.0],
            )
            .unwrap(),
        )
        .unwrap();
    let gas = ClosedGasParcel::new(pure_b, point(400.0, 100_000.0), 0.001).unwrap();
    let flux = law.sample_to_vacuum(&gas, AREA, 900.0, || false).unwrap();
    assert_eq!(flux.transmitted.species.len(), 1);
    assert_eq!(flux.wall.reflected_species.len(), 1);
    assert_eq!(
        flux.transmitted.species[0].species,
        SpeciesId::new("B").unwrap()
    );
    assert_eq!(
        flux.wall.reflected_species[0].species,
        flux.transmitted.species[0].species
    );
    near(
        flux.transmitted.molar_rate_mol_s,
        (18.0 / 35.0) * incident_moles(100_000.0, 400.0, 0.004),
    );
}
