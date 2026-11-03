//! Independent Maxwellian moments; authored gases are equation tests, not data.
use super::*;
use crate::thermochemical_gas::{
    ThermochemicalGasModel,
    energy::parcel::transfer::species::{
        effusion::{EffusiveRegime, EffusiveSpeciesTransmission, EffusiveSupport},
        network::{
            SpeciesNetworkSolve,
            transport::{
                GasTransportPaths, GasTransportSolve,
                effusion::vacuum::{VacuumEffusionLink, step_candidate},
            },
        },
    },
    energy::{GasEnergySolve, parcel::ClosedGasParcel},
    tests::{point, species_with_energy_reference, transport},
};
use fs_thermochem::{
    Composition, CompositionBasis, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
    mixture::FrozenIdealGasMixtureModelV1,
};
use std::collections::BTreeMap;

const AREA: f64 = 1e-8;
fn gas(datum: f64, temperature: f64, pressure: f64) -> ClosedGasParcel {
    let mixture = FrozenIdealGasMixtureModelV1::new(
        vec![
            species_with_energy_reference("A", 0.004, 2.5, 0.0, datum),
            species_with_energy_reference("B", 0.028, 3.5, 0.001, datum),
        ],
        Composition::new(CompositionBasis::MoleFraction, vec![0.25, 0.75]).unwrap(),
    )
    .unwrap();
    ClosedGasParcel::new(
        ThermochemicalGasModel::new(
            mixture,
            vec![transport("A", 1.7e-5, 0.025), transport("B", 1.7e-5, 0.025)],
            None,
        )
        .unwrap(),
        point(temperature, pressure),
        0.001,
    )
    .unwrap()
}
fn qualified(gas: &ClosedGasParcel, tau: [f64; 2]) -> EffusiveMomentumInterface {
    let rows = gas
        .model()
        .thermodynamics()
        .components()
        .iter()
        .zip(tau)
        .map(|(s, probability)| {
            EffusiveSpeciesTransmission::try_new(
                s.species().clone(),
                probability,
                "synthetic reciprocal energy-independent transmission",
            )
            .unwrap()
        })
        .collect();
    let interface = EffusiveInterface::try_new(
        gas.model().clone(),
        rows,
        EffusiveSupport {
            temperature_k: (200.0, 2000.0),
            pressure_pa: (1000.0, 1e7),
        },
        EffusiveRegime::try_new(
            1e-9,
            1e-12,
            1e-7,
            "synthetic thin ballistic aperture for the complete gas support",
        )
        .unwrap(),
    )
    .unwrap();
    EffusiveMomentumInterface::try_new(
        interface,
        "synthetic angle-independent nondeflecting transmission for every species and state",
    )
    .unwrap()
}
fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 8e-12 * actual.abs().max(expected.abs()).max(1e-250),
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn species_normal_moments_recover_half_partial_pressure_without_an_energy_speed() {
    for datum in [0.0, 1e12] {
        let parcel = gas(datum, 400.0, 100_000.0);
        let qualification = qualified(&parcel, [0.8, 0.3]);
        let original = parcel
            .effusive_flux_to_vacuum(qualification.interface(), AREA, || false)
            .unwrap();
        let moments = qualification.vacuum_flux(&original, || false).unwrap();
        assert!(std::ptr::eq(moments.flux(), &original));
        assert!(std::ptr::eq(moments.qualification(), &qualification));
        for (index, (molar_mass, fraction, tau)) in [(0.004, 0.25, 0.8), (0.028, 0.75, 0.3)]
            .into_iter()
            .enumerate()
        {
            let row = &moments.species()[index];
            let expected_speed = (core::f64::consts::PI * R * 400.0 / (2.0 * molar_mass)).sqrt();
            near(row.outward.mean_normal_speed_m_s().unwrap(), expected_speed);
            near(
                row.outward.normal_momentum_rate_n(),
                AREA * tau * fraction * 100_000.0 / 2.0,
            );
            assert!(std::ptr::eq(
                row.outward.beam(),
                &original.species[index].outward
            ));
            near(
                row.outward.impulse_for_mass_kg(1e-9).unwrap(),
                1e-9 * expected_speed,
            );
        }
        near(
            moments.normal_momentum_rate_n(),
            AREA * 100_000.0 * (0.8 * 0.25 + 0.3 * 0.75) / 2.0,
        );
    }
}

#[test]
fn zero_net_mass_still_transmits_pressure_stress_in_both_directions() {
    let left = gas(0.0, 400.0, 100_000.0);
    let qualification = qualified(&left, [1.0, 1.0]);
    for (temperature, pressure) in [(400.0, 100_000.0), (1600.0, 200_000.0)] {
        let right = gas(0.0, temperature, pressure);
        let original = left
            .effusive_flux_to(&right, qualification.interface(), AREA, || false)
            .unwrap();
        let moments = qualification.reciprocal_flux(&original, || false).unwrap();
        assert!(original.net_mass_rate_left_to_right_kg_s.abs() < 1e-20);
        near(
            moments.normal_momentum_stress_n(),
            AREA * (100_000.0 + pressure) / 2.0,
        );
        for (row, fraction) in moments.species().iter().zip([0.25, 0.75]) {
            near(
                row.left_to_right.normal_momentum_rate_n(),
                AREA * fraction * 100_000.0 / 2.0,
            );
            near(
                row.right_to_left.normal_momentum_rate_n(),
                AREA * fraction * pressure / 2.0,
            );
        }
        let reverse = right
            .effusive_flux_to(&left, qualification.interface(), AREA, || false)
            .unwrap();
        near(
            qualification
                .reciprocal_flux(&reverse, || false)
                .unwrap()
                .normal_momentum_stress_n(),
            moments.normal_momentum_stress_n(),
        );
    }
}

#[test]
fn finite_momentum_uses_each_original_accepted_species_withdrawal() {
    let parcel = gas(1e12, 400.0, 100_000.0);
    let qualification = qualified(&parcel, [0.8, 0.3]);
    let parcels = BTreeMap::from([(0, parcel)]);
    let endpoints = BTreeMap::from([(
        0,
        GasEnergySolve {
            minimum_temperature_k: 200.0,
            maximum_temperature_k: 1800.0,
            energy_tolerance_j_kg: 1e-7,
            maximum_iterations: 128,
        },
    )]);
    let candidate = step_candidate(
        &parcels,
        GasTransportPaths::default(),
        &[VacuumEffusionLink {
            id: 7,
            parcel: 0,
            interface: qualification.interface(),
            area_m2: AREA,
        }],
        0.01,
        &endpoints,
        GasTransportSolve {
            exchange: SpeciesNetworkSolve {
                maximum_species_fraction: 0.1,
                mass_tolerance_kg: 1e-13,
                energy_tolerance_j: 1e-8,
            },
            entropy_tolerance_j_k: 1e-10,
        },
        || false,
    )
    .unwrap();
    let original = &candidate.step().outlets[0].flux;
    let moments = qualification.vacuum_flux(original, || false).unwrap();
    for (row, (molar_mass, fraction, tau)) in moments
        .species()
        .iter()
        .zip([(0.004, 0.25, 0.8), (0.028, 0.75, 0.3)])
    {
        let account = candidate
            .step()
            .exchange
            .boundaries
            .iter()
            .find(|account| &account.transfer.species == row.species)
            .unwrap();
        assert_eq!(account.id, 7);
        let mass = account.transfer.outgoing_mass_kg;
        near(
            mass,
            AREA * tau
                * fraction
                * 100_000.0
                * (molar_mass / (2.0 * core::f64::consts::PI * R * 400.0)).sqrt()
                * 0.01,
        );
        near(
            row.outward.impulse_for_mass_kg(mass).unwrap(),
            mass * (core::f64::consts::PI * R * 400.0 / (2.0 * molar_mass)).sqrt(),
        );
        // The extensive impulse follows the admitted mass even if a consumer
        // splits that original parcel into two accepted accounting windows.
        near(
            row.outward.impulse_for_mass_kg(mass / 2.0).unwrap() * 2.0,
            row.outward.impulse_for_mass_kg(mass).unwrap(),
        );
    }
}

#[test]
fn closed_paths_have_no_invented_velocity_or_impulse() {
    let parcel = gas(0.0, 400.0, 100_000.0);
    let qualification = qualified(&parcel, [0.0, 0.0]);
    for area in [0.0, AREA] {
        let original = parcel
            .effusive_flux_to_vacuum(qualification.interface(), area, || false)
            .unwrap();
        let moments = qualification.vacuum_flux(&original, || false).unwrap();
        assert_eq!(moments.normal_momentum_rate_n(), 0.0);
        for row in moments.species() {
            assert_eq!(row.outward.mean_normal_speed_m_s(), None);
            assert_eq!(row.outward.impulse_for_mass_kg(0.0).unwrap(), 0.0);
            assert!(row.outward.impulse_for_mass_kg(f64::MIN_POSITIVE).is_err());
        }
    }
}

#[test]
fn angular_qualification_stays_bound_and_cancellation_is_atomic() {
    let parcel = gas(0.0, 400.0, 100_000.0);
    let qualification = qualified(&parcel, [1.0, 1.0]);
    let other = qualified(&parcel, [1.0, 1.0]);
    let original = parcel
        .effusive_flux_to_vacuum(qualification.interface(), AREA, || false)
        .unwrap();
    assert!(other.vacuum_flux(&original, || false).is_err());
    assert!(
        qualification
            .clone()
            .vacuum_flux(&original, || false)
            .is_ok()
    );
    assert!(qualification.vacuum_flux(&original, || true).is_err());
    let mut changed_species = original.clone();
    changed_species.species[0].species = SpeciesId::new("unrelated").unwrap();
    assert!(
        qualification
            .vacuum_flux(&changed_species, || false)
            .is_err()
    );
    let mut duplicate = original.clone();
    duplicate.species.push(duplicate.species[0].clone());
    assert!(qualification.vacuum_flux(&duplicate, || false).is_err());
    assert!(EffusiveMomentumInterface::try_new(qualification.interface().clone(), " ").is_err());
    assert!(
        EffusiveMomentumInterface::try_new(qualification.interface().clone(), "x".repeat(4097))
            .is_err()
    );
    let moments = qualification.vacuum_flux(&original, || false).unwrap();
    for mass in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            moments.species()[0]
                .outward
                .impulse_for_mass_kg(mass)
                .is_err()
        );
    }
}

#[test]
fn normal_moment_keeps_representable_results_across_extreme_binary_scales() {
    // Isolate arithmetic with admitted moment values: pi*1e308 overflows,
    // although both the square-root speed and its product with mass are finite.
    let beam = EffusiveBeam {
        mass_rate_kg_s: 1e-150,
        molar_rate_mol_s: 1e-150,
        additional_specific_energy_j_kg: Some(-1e308),
        ..EffusiveBeam::default()
    };
    let moment = momentum(&beam).unwrap();
    near(
        moment.mean_normal_speed_m_s().unwrap(),
        core::f64::consts::PI.sqrt() * 1e154,
    );
    near(
        moment.normal_momentum_rate_n(),
        core::f64::consts::PI.sqrt() * 1e4,
    );
    near(
        moment.impulse_for_mass_kg(1e-310).unwrap(),
        1e-310 * (core::f64::consts::PI.sqrt() * 1e154),
    );
    let overflow = EffusiveBeam {
        mass_rate_kg_s: 1e308,
        ..beam
    };
    assert!(momentum(&overflow).is_err());
    let tiny = EffusiveBeam {
        mass_rate_kg_s: 1e300,
        additional_specific_energy_j_kg: Some(-f64::from_bits(1)),
        ..beam
    };
    near(
        momentum(&tiny).unwrap().normal_momentum_rate_n(),
        (core::f64::consts::PI.sqrt() * f64::from_bits(1).sqrt()) * 1e300,
    );
}
