//! Manufactured source/branch oracles for property-only real-fluid evaluations.
use super::super::tests::{ternary_bindings, ternary_source};
use super::*;
use crate::thermochemical_gas::tests::{model, point};
use fs_thermochem::{Composition, CompositionBasis};

fn interactions() -> Vec<f64> {
    vec![0.0, 0.02, 0.12, 0.02, 0.0, -0.03, 0.12, -0.03, 0.0]
}

#[test]
fn extra_species_bank_retains_actual_mass_composition_and_both_pair_axes() {
    let source = ternary_source();
    let bindings = ternary_bindings();
    let bank =
        ThermochemicalPengRobinsonModel::new(source.clone(), bindings.clone(), interactions())
            .unwrap();
    // Removing the middle species exercises both non-contiguous matrix axes.
    let active = ThermochemicalPengRobinsonModel::new(
        source
            .with_composition(
                Composition::new(CompositionBasis::MassFraction, vec![0.7, 0.0, 0.3]).unwrap(),
            )
            .unwrap(),
        vec![bindings[0].clone(), bindings[2].clone()],
        vec![0.0, 0.12, 0.12, 0.0],
    )
    .unwrap();
    let context = point(350.0, 5.0e6).with("exposure", 7.5).unwrap();
    let bulk = active
        .evaluate_real_gas(&context, PengRobinsonPhaseV1::Vapor)
        .unwrap();
    let saved = bulk.clone();
    // Common-bank trials and a larger source bank must resolve identical states.
    let direct = active
        .evaluate_temperature_from_snapshot(
            &bulk,
            Temperature::new(375.0),
            PengRobinsonPhaseV1::Vapor,
        )
        .unwrap();
    let film = bank
        .evaluate_temperature_from_snapshot(
            &bulk,
            Temperature::new(375.0),
            PengRobinsonPhaseV1::Vapor,
        )
        .unwrap();
    assert_eq!(film, direct);
    assert!(std::ptr::eq(
        film.constitutive_law(),
        bulk.constitutive_law()
    ));
    assert_eq!(film.temperature_k(), 375.0);
    assert_eq!(film.pressure_pa(), bulk.pressure_pa());
    assert_eq!(film.source().components().len(), 2);
    let receipt = film.source().thermodynamics().receipt();
    let original = bulk.source().thermodynamics().receipt();
    assert_eq!(receipt.declared_basis(), CompositionBasis::MassFraction);
    assert_eq!(
        receipt.declared_fraction_bits(),
        original.declared_fraction_bits()
    );
    assert_eq!(receipt.mole_fraction_bits(), original.mole_fraction_bits());
    assert_eq!(receipt.species(), original.species());
    assert_eq!(film.source().source_packs(), bulk.source().source_packs());
    let expected_context = context.with("temperature", 375.0).unwrap();
    for (property, original) in film
        .source()
        .components()
        .iter()
        .zip(bulk.source().components())
    {
        assert_eq!(property.query(), expected_context);
        assert_eq!(property.card_identity(), original.card_identity());
        for (a, b) in property.properties().iter().zip(original.properties()) {
            assert_eq!(a.answer().receipt.selected, b.answer().receipt.selected);
        }
    }
    let expected = active
        .evaluate_real_gas(&expected_context, PengRobinsonPhaseV1::Vapor)
        .unwrap();
    assert_eq!(film.root(), expected.root());
    assert_eq!(
        film.transport_state().unwrap(),
        expected.transport_state().unwrap()
    );
    assert_eq!(
        film.volumetric_expansion_per_k(),
        expected.volumetric_expansion_per_k()
    );
    assert_ne!(film.density_kg_m3(), film.source().state().density);

    // Inactive-only law differences must not reactivate B or perturb the film.
    let mut unused_pairs = interactions();
    unused_pairs[1] = 0.22;
    unused_pairs[3] = 0.22;
    unused_pairs[5] = 0.17;
    unused_pairs[7] = 0.17;
    let mut unused_critical = bindings.clone();
    unused_critical[1].critical_temperature = Temperature::new(180.0);
    let inactive_changes =
        ThermochemicalPengRobinsonModel::new(source.clone(), unused_critical, unused_pairs)
            .unwrap();
    assert_eq!(
        inactive_changes
            .evaluate_temperature_from_snapshot(
                &bulk,
                Temperature::new(375.0),
                PengRobinsonPhaseV1::Vapor,
            )
            .unwrap(),
        film
    );

    let mut changed_pairs = interactions();
    changed_pairs[2] = 0.13;
    changed_pairs[6] = 0.13;
    let changed_pair =
        ThermochemicalPengRobinsonModel::new(source.clone(), bindings.clone(), changed_pairs)
            .unwrap();
    let mut changed_critical = bindings;
    changed_critical[2].critical_temperature = Temperature::new(305.0);
    let changed_law =
        ThermochemicalPengRobinsonModel::new(source, changed_critical, interactions()).unwrap();
    for changed in [&changed_pair, &changed_law] {
        assert!(
            changed
                .evaluate_temperature_from_snapshot(
                    &bulk,
                    Temperature::new(375.0),
                    PengRobinsonPhaseV1::Vapor,
                )
                .is_err()
        );
    }
    assert_eq!(bulk, saved);
}

#[test]
fn explicit_branch_retains_liquid_or_vapor_and_refuses_silent_reselection() {
    let source = model(vec![0.75, 0.25], 0.0)
        .with_composition(Composition::new(CompositionBasis::MoleFraction, vec![1.0, 0.0]).unwrap())
        .unwrap();
    // Twice the methane critical temperature maps the original owner's
    // subcritical 150 K / 1 MPa three-root oracle to source-supported 300 K.
    let bound = ThermochemicalPengRobinsonModel::new(
        source,
        vec![PengRobinsonSpeciesBinding {
            species: SpeciesId::new("A").unwrap(),
            critical_temperature: Temperature::new(381.128),
            critical_pressure: Pressure::new(4.5992e6),
            acentric_factor: 0.01142,
        }],
        vec![],
    )
    .unwrap();
    for (phase, wrong) in [
        (PengRobinsonPhaseV1::Vapor, PengRobinsonPhaseV1::Liquid),
        (PengRobinsonPhaseV1::Liquid, PengRobinsonPhaseV1::Vapor),
    ] {
        let bulk = bound
            .evaluate_real_gas(&point(300.0, 1.0e6), phase)
            .unwrap();
        assert_eq!(bulk.mechanics().phase_roots().len(), 3);
        let film = bound
            .evaluate_temperature_from_snapshot(&bulk, Temperature::new(302.0), phase)
            .unwrap();
        assert_eq!(film.mechanics().phase_roots().len(), 3);
        assert_eq!(film.root(), film.mechanics().select(phase));
        assert_ne!(film.root(), film.mechanics().select(wrong));
        assert_eq!(film.pressure_pa(), bulk.pressure_pa());
        for refused in [wrong, PengRobinsonPhaseV1::Stable] {
            assert!(
                bound
                    .evaluate_temperature_from_snapshot(&bulk, Temperature::new(302.0), refused,)
                    .is_err()
            );
        }
    }
}

#[test]
fn shared_mechanics_still_admit_original_sources_and_property_domains() {
    let source = ternary_source();
    let bank =
        ThermochemicalPengRobinsonModel::new(source, ternary_bindings(), interactions()).unwrap();
    let bulk = bank
        .evaluate_real_gas(&point(350.0, 5e6), PengRobinsonPhaseV1::Vapor)
        .unwrap();
    let saved = bulk.clone();
    // Sharing mechanical storage cannot exempt the caloric/transport sources.
    let mut changed_source = bank.clone();
    changed_source.source.transport[0] =
        crate::thermochemical_gas::tests::transport("A", 2.0e-5, 0.03);
    assert!(Arc::ptr_eq(&changed_source.mechanics, &bulk.law));
    assert!(
        changed_source
            .evaluate_temperature_from_snapshot(
                &bulk,
                Temperature::new(375.0),
                PengRobinsonPhaseV1::Vapor,
            )
            .is_err()
    );
    for temperature in [f64::NAN, f64::INFINITY, 0.0, -1.0, 199.0, 2001.0] {
        assert!(
            bank.evaluate_temperature_from_snapshot(
                &bulk,
                Temperature::new(temperature),
                PengRobinsonPhaseV1::Vapor,
            )
            .is_err()
        );
    }
    assert!(
        bank.evaluate_temperature_from_snapshot(
            &bulk,
            Temperature::new(350.0),
            PengRobinsonPhaseV1::Stable,
        )
        .is_err()
    );
    assert_eq!(bulk, saved);
}
