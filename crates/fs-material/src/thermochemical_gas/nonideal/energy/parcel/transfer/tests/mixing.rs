//! Independent mixture-energy oracle and source/transaction boundaries.

use super::*;
use crate::fluid_transport::FluidTransport;
use crate::thermochemical_gas::tests::{
    species_in_domain, species_with_energy_reference, transport,
};
use fs_qty::Amount;
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::{Composition, CompositionBasis};

fn source(x_a: f64, references: [f64; 2]) -> ThermochemicalGasModel {
    ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![
                species_with_energy_reference("A", 0.028, 3.5, 0.001, references[0]),
                species_with_energy_reference("B", 0.004, 2.5, 0.0, references[1]),
            ],
            Composition::new(CompositionBasis::MoleFraction, vec![x_a, 1.0 - x_a]).unwrap(),
        )
        .unwrap(),
        vec![transport("A", 1.7e-5, 0.025), transport("B", 2.0e-5, 0.15)],
        None,
    )
    .unwrap()
}

fn bindings(source: &ThermochemicalGasModel) -> Vec<PengRobinsonSpeciesBinding> {
    let original = bound();
    source
        .thermodynamics()
        .components()
        .iter()
        .map(|species| {
            let i = original
                .source()
                .thermodynamics()
                .components()
                .iter()
                .position(|candidate| candidate.species() == species.species())
                .unwrap();
            let data = original.mechanics().components()[i];
            PengRobinsonSpeciesBinding {
                species: species.species().clone(),
                critical_temperature: data.critical_temperature,
                critical_pressure: data.critical_pressure,
                acentric_factor: data.acentric_factor,
            }
        })
        .collect()
}

fn real_model(source: ThermochemicalGasModel, interaction: f64) -> ThermochemicalPengRobinsonModel {
    let data = bindings(&source);
    // Empty pure catalogs know their own diagonal only, never the absent AB pair.
    let interactions = if data.len() == 1 {
        vec![]
    } else {
        vec![0.0, interaction, interaction, 0.0]
    };
    ThermochemicalPengRobinsonModel::new(source, data, interactions).unwrap()
}

fn at(
    model: ThermochemicalPengRobinsonModel,
    temperature: f64,
    pressure: f64,
    owner: f64,
) -> IsochoricRealGasParcel {
    IsochoricRealGasParcel::new(
        model,
        point(temperature, pressure).with("owner", owner).unwrap(),
        0.001,
        PengRobinsonPhaseV1::Vapor,
    )
    .unwrap()
}

fn pair(references: [f64; 2]) -> [IsochoricRealGasParcel; 2] {
    [
        at(real_model(source(0.75, references), 0.03), 300.0, 5e6, 11.0),
        at(
            real_model(source(0.35, references), 0.03),
            350.0,
            3.8e6,
            22.0,
        ),
    ]
}

// Independent extensive conversion for the manufactured two-species fixture.
fn species_mass(gas: &IsochoricRealGasParcel, name: &str) -> f64 {
    let id = SpeciesId::new(name).unwrap();
    gas.model()
        .source()
        .thermodynamics()
        .components()
        .iter()
        .zip(gas.model().mechanics().mole_fractions())
        .find(|(species, _)| species.species() == &id)
        .map_or(0.0, |(species, x)| {
            gas.mass_kg() * x * species.molar_mass().value()
                / gas.model().mechanics().molar_mass().value()
        })
}

#[test]
fn different_compositions_match_independent_pr_mixing_energy_and_continue() {
    let [donor, receiver] = pair([0.0; 2]);
    let candidate = donor
        .transfer_mass_to(&receiver, 0.001, controls(), || false)
        .unwrap();
    let step = candidate.step();
    let after = [candidate.donor(), candidate.receiver()];
    // Independent 80-digit Decimal solution of the original PR/NASA laws,
    // conserving each species and U_d+U_r with dm*h_d transported once.
    close(after[0].mass_kg(), 0.046_889_003_167_633_395, 1e-14);
    close(after[0].state().temperature_k(), 297.148_041_580_514, 2e-8);
    close(after[0].state().pressure_pa(), 4_845_402.103_709_031, 5e-4);
    close(after[1].mass_kg(), 0.017_555_285_197_992_412, 1e-14);
    close(
        after[1].state().temperature_k(),
        352.535_331_755_326_26,
        2e-8,
    );
    close(
        after[1].state().pressure_pa(),
        3_956_630.450_099_272_6,
        5e-4,
    );
    close(
        after[1].model().mechanics().mole_fractions()[0],
        0.363_169_903_527_094_3,
        1e-14,
    );
    close(step.transferred_enthalpy_j, 346.804_450_075_851_95, 1e-8);
    close(step.donor_energy_change_j, -346.804_450_075_851_95, 1e-7);
    close(step.receiver_energy_change_j, 346.804_450_075_851_95, 1e-7);
    let initial_u = donor.mass_kg() * donor.state().specific_internal_energy_j_kg()
        + receiver.mass_kg() * receiver.state().specific_internal_energy_j_kg();
    let final_u = after
        .iter()
        .map(|gas| gas.mass_kg() * gas.state().specific_internal_energy_j_kg())
        .sum::<f64>();
    close(initial_u, 18_682.263_985_181_478, 1e-7);
    close(final_u, initial_u, controls().energy_tolerance_j);
    assert_eq!(step.species.len(), 2);
    for (i, name) in ["A", "B"].into_iter().enumerate() {
        assert_eq!(step.species[i].species, SpeciesId::new(name).unwrap());
        let expected = if i == 0 {
            0.000_954_545_454_545_454_6
        } else {
            0.000_045_454_545_454_545_455
        };
        close(step.species[i].mass_kg, expected, 1e-17);
        close(
            species_mass(after[0], name) - species_mass(&donor, name),
            -expected,
            1e-14,
        );
        close(
            species_mass(after[1], name) - species_mass(&receiver, name),
            expected,
            1e-14,
        );
    }
    close(
        species_mass(after[1], "A"),
        0.014_038_561_175_539_458,
        1e-14,
    );
    close(
        species_mass(after[1], "B"),
        0.003_516_724_022_452_954_4,
        1e-14,
    );
    assert!(step.maximum_species_mass_residual_kg <= controls().mass_tolerance_kg);
    assert!(step.energy_residual_j.abs() <= controls().energy_tolerance_j);
    for (gas, owner) in after.into_iter().zip([11.0, 22.0]) {
        assert_eq!(gas.volume_m3(), 0.001);
        assert_eq!(gas.density_kg_m3(), gas.mass_kg() / gas.volume_m3());
        let state = gas.state();
        assert_eq!(state.density_kg_m3(), gas.density_kg_m3());
        assert_eq!(
            state.transport_state().unwrap().density_kg_m3(),
            gas.density_kg_m3()
        );
        assert_eq!(
            state.transport_state().unwrap().specific_heat_j_kg_k(),
            state.specific_heat_cp_j_kg_k()
        );
        assert!(state.specific_heat_cv_j_kg_k() > 0.0 && state.sound_speed_m_s() > 0.0);
        assert!((state.source().state().density - state.density_kg_m3()).abs() > 1e-3);
        for component in state.source().components() {
            assert_eq!(component.query().axes().get("owner"), Some(&owner));
        }
    }
    // The returned complete receiver immediately participates in the original
    // heat/work owner with its new composition and its original source cards.
    let moved = candidate
        .receiver()
        .midpoint_heat_volume(0.001_001, 0.5, energy(), || false)
        .unwrap();
    assert_eq!(moved.endpoint().model(), candidate.receiver().model());
    assert_eq!(moved.endpoint().mass_kg(), candidate.receiver().mass_kg());
    close(
        moved.step().stored_energy_change_j + moved.step().work_by_gas_j,
        0.5,
        1e-7,
    );
}

#[test]
fn distinct_huge_species_formation_datums_do_not_change_mixing_physics() {
    let [donor, receiver] = pair([0.0; 2]);
    let expected = donor
        .transfer_mass_to(&receiver, 0.001, controls(), || false)
        .unwrap();
    for references in [[1e24, -2e24], [-3e24, 4e24]] {
        let [donor, receiver] = pair(references);
        let actual = donor
            .transfer_mass_to(&receiver, 0.001, controls(), || false)
            .unwrap();
        for (actual, expected) in [actual.donor(), actual.receiver()]
            .into_iter()
            .zip([expected.donor(), expected.receiver()])
        {
            assert_eq!(
                actual.state().temperature_k(),
                expected.state().temperature_k()
            );
            assert_eq!(actual.state().pressure_pa(), expected.state().pressure_pa());
            assert_eq!(actual.density_kg_m3(), expected.density_kg_m3());
            assert_eq!(actual.model().mechanics(), expected.model().mechanics());
        }
        assert_eq!(actual.step().species, expected.step().species);
        assert_eq!(
            actual.step().energy_residual_j,
            expected.step().energy_residual_j
        );
        // Absolute donor u has lost the entire physical cooling signal; the
        // native increment still resolves it without subtracting these values.
        assert_eq!(
            actual.donor().state().specific_internal_energy_j_kg(),
            donor.state().specific_internal_energy_j_kg()
        );
        assert!(actual.donor().state().temperature_k() < donor.state().temperature_k() - 2.0);
    }
}

#[test]
fn incoming_species_keeps_known_pair_data_and_missing_cross_pair_refuses() {
    let base = source(0.75, [0.0; 2]);
    let pure_a = base
        .with_mole_amounts(&[Amount::new(1.0), Amount::new(0.0)])
        .unwrap();
    let pure_b = base
        .with_mole_amounts(&[Amount::new(0.0), Amount::new(1.0)])
        .unwrap();
    let donor = at(real_model(base, 0.03), 300.0, 5e6, 11.0);
    let receiver = at(real_model(pure_b.clone(), 0.03), 350.0, 3.8e6, 22.0);
    let candidate = donor
        .transfer_mass_to(&receiver, 0.001, controls(), || false)
        .unwrap();
    assert_eq!(receiver.model().mechanics().components().len(), 1);
    assert_eq!(
        candidate.receiver().model().mechanics().components().len(),
        2
    );
    assert_eq!(
        candidate
            .receiver()
            .model()
            .mechanics()
            .binary_interactions(),
        &[0.0, 0.03, 0.03, 0.0]
    );
    close(
        species_mass(candidate.receiver(), "A"),
        candidate.step().species[0].mass_kg,
        1e-14,
    );
    // A second transfer can now use the combined original source inventory.
    let reversed = candidate
        .receiver()
        .transfer_mass_to(candidate.donor(), 0.0001, controls(), || false)
        .unwrap();
    assert!(reversed.step().maximum_species_mass_residual_kg <= controls().mass_tolerance_kg);
    assert!(reversed.step().energy_residual_j.abs() <= controls().energy_tolerance_j);

    let donor = at(real_model(pure_a, 0.03), 300.0, 5e6, 11.0);
    let saved = [donor.clone(), receiver.clone()];
    let zero = donor
        .transfer_mass_to(&receiver, 0.0, controls(), || false)
        .unwrap();
    assert_eq!([zero.donor(), zero.receiver()], [&donor, &receiver]);
    assert!(zero.step().species.is_empty());
    assert!(matches!(
        donor.transfer_mass_to(&receiver, 0.001, controls(), || false),
        Err(GasEnergyError::Invalid(
            "mixed gas needs explicit Peng-Robinson interaction data for every new species pair"
        ))
    ));
    assert_eq!([donor, receiver], saved);
}

#[test]
fn shared_critical_interaction_and_original_source_conflicts_refuse_atomically() {
    let [donor, receiver] = pair([0.0; 2]);
    let mut critical = bindings(receiver.model().source());
    critical[0].critical_temperature = Temperature::new(191.0);
    let wrong_critical = ThermochemicalPengRobinsonModel::new(
        receiver.model().source().clone(),
        critical,
        vec![0.0, 0.03, 0.03, 0.0],
    )
    .unwrap();
    let wrong_transport = ThermochemicalGasModel::new(
        receiver.model().source().thermodynamics().clone(),
        vec![transport("A", 1.8e-5, 0.025), transport("B", 2.0e-5, 0.15)],
        None,
    )
    .unwrap();
    let incompatible = [
        wrong_critical,
        real_model(receiver.model().source().clone(), 0.04),
        real_model(model(vec![0.35, 0.65], 0.002), 0.03),
        real_model(wrong_transport, 0.03),
    ];
    let conflicts = [
        "shared gas species has conflicting Peng-Robinson critical properties",
        "shared gas species pair has conflicting Peng-Robinson interaction data",
        "shared gas species has conflicting thermochemical or transport source laws",
        "shared gas species has conflicting thermochemical or transport source laws",
    ];
    for (wrong, reason) in incompatible.into_iter().zip(conflicts) {
        let receiver = at(wrong, 350.0, 3.8e6, 22.0);
        let saved = [donor.clone(), receiver.clone()];
        assert!(matches!(
            donor.transfer_mass_to(&receiver, 0.001, controls(), || false),
            Err(GasEnergyError::Invalid(actual)) if actual == reason
        ));
        assert_eq!([&donor, &receiver], [&saved[0], &saved[1]]);
    }
}

#[test]
fn incoming_source_must_cover_receiving_anchor_and_late_refusals_are_atomic() {
    let [donor, receiver] = pair([0.0; 2]);
    let saved = [donor.clone(), receiver.clone()];
    let mut polls = 0;
    let accepted = donor
        .transfer_mass_to(&receiver, 0.001, controls(), || {
            polls += 1;
            false
        })
        .unwrap();
    let mut seen = 0;
    assert!(matches!(
        donor.transfer_mass_to(&receiver, 0.001, controls(), || {
            seen += 1;
            seen == polls
        }),
        Err(GasEnergyError::Cancelled)
    ));
    assert_eq!(seen, polls);
    let mut exhausted = controls();
    exhausted.receiver.maximum_iterations = 1;
    assert!(matches!(
        donor.transfer_mass_to(&receiver, 0.001, exhausted, || false),
        Err(GasEnergyError::NotConverged { .. })
    ));
    assert_eq!([&donor, &receiver], [&saved[0], &saved[1]]);
    assert_eq!(
        donor
            .transfer_mass_to(&receiver, 0.001, controls(), || false)
            .unwrap(),
        accepted
    );

    let incoming = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![
                species_in_domain("A", 0.028, 3.5, 0.001, 0.0, [200.0, 320.0]),
                species_with_energy_reference("B", 0.004, 2.5, 0.0, 0.0),
            ],
            Composition::new(CompositionBasis::MoleFraction, vec![0.75, 0.25]).unwrap(),
        )
        .unwrap(),
        vec![transport("A", 1.7e-5, 0.025), transport("B", 2.0e-5, 0.15)],
        None,
    )
    .unwrap();
    let donor = at(real_model(incoming, 0.03), 300.0, 5e6, 11.0);
    let pure_b = source(0.35, [0.0; 2])
        .with_mole_amounts(&[Amount::new(0.0), Amount::new(1.0)])
        .unwrap();
    let receiver = at(real_model(pure_b, 0.03), 350.0, 3.8e6, 22.0);
    let mut domain_controls = controls();
    domain_controls.donor.maximum_temperature_k = 310.0;
    domain_controls.receiver.minimum_temperature_k = 330.0;
    domain_controls.receiver.maximum_temperature_k = 400.0;
    let saved = [donor.clone(), receiver.clone()];
    assert!(matches!(
        donor.transfer_mass_to(&receiver, 0.001, domain_controls, || false),
        Err(GasEnergyError::Material(_))
    ));
    assert_eq!([donor, receiver], saved);
}
