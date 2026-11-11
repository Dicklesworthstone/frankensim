//! Original sources, frozen populations and actual temperatures survive refresh.
use super::*;
use fs_material::ionized_gas::{CarrierMomentumSource, IonizedGasRefreshError};
use fs_qty::{Temperature, Volume};
use fs_thermochem::equilibrium::composition::state::FrozenIonizedGasStateV1;

fn heated(old: &FrozenIonizedGasStateV1, th: f64, te: f64) -> FrozenIonizedGasStateV1 {
    FrozenIonizedGasStateV1::new(
        old.thermodynamics(),
        old.amounts(),
        old.volume(),
        Temperature::new(th),
        Temperature::new(te),
        1e-12,
        || false,
    )
    .unwrap()
}
fn source_with_card(
    state: &FrozenIonizedGasStateV1,
) -> (IonizedGasConductivity, fs_matdb::MaterialCard) {
    let (ion, card) = collision_source_and_card(
        "H+",
        &[
            ("T_h", state.heavy_temperature().value()),
            ("T_e", state.electron_temperature().value()),
        ],
        state.pressure().value(),
        1e9,
    );
    let electron = electron_relaxation(state, 2e-20);
    (
        IonizedGasConductivity::from_state(state.clone(), vec![ion, electron], limits(), || false)
            .unwrap(),
        card,
    )
}

#[test]
fn original_cross_sections_are_reintegrated_at_actual_temperature_and_density() {
    let state = frozen_state(0.001, 6000.0, 12000.0, 1.0);
    let (original, card) = source_with_card(&state);
    let next = heated(&state, 6500.0, 15000.0);
    let refreshed = original.at_frozen_state(&next, &[&card], || false).unwrap();
    assert_eq!(refreshed.state(), &next);
    assert_eq!(refreshed.limits(), original.limits());
    assert!(refreshed.equilibrium().is_none());
    assert_ne!(refreshed.conductivity_s_m(), original.conductivity_s_m());
    for (before, after) in original.carriers().iter().zip(refreshed.carriers()) {
        assert_eq!(before.species(), after.species());
        assert_eq!(before.mass_kg(), after.mass_kg());
        assert_eq!(
            before.number_density_per_m3(),
            after.number_density_per_m3()
        );
        relative(
            after.relaxation().pressure_pa(),
            next.pressure().value(),
            2e-15,
        );
        match (before.relaxation().source(), after.relaxation().source()) {
            (
                CarrierMomentumSource::ElectronNeutral {
                    contributions: old, ..
                },
                CarrierMomentumSource::ElectronNeutral {
                    state: bound,
                    contributions: new,
                    ..
                },
            ) => {
                assert_eq!(bound, &next);
                assert_eq!(old[0].evaluation().source(), new[0].evaluation().source());
                assert_eq!(
                    old[0].evaluation().controls(),
                    new[0].evaluation().controls()
                );
                assert_eq!(
                    new[0].evaluation().electron_temperature(),
                    next.electron_temperature()
                );
                relative(
                    after.relaxation().rate_per_s() / before.relaxation().rate_per_s(),
                    (15000.0_f64 / 12000.0).sqrt(),
                    2e-8,
                );
            }
            (
                CarrierMomentumSource::Material {
                    material: old,
                    property: key,
                },
                CarrierMomentumSource::Material {
                    material: new,
                    property: refreshed_key,
                },
            ) => {
                assert_eq!(old.card_identity(), new.card_identity());
                assert_eq!(key, refreshed_key);
                assert_eq!(new.state_coordinate("T_h"), Some(6500.0));
                assert_eq!(new.state_coordinate("T_e"), Some(15000.0));
                assert_eq!(new.axis_quantities(), old.axis_quantities());
            }
            _ => panic!("refresh changed the original momentum-source kind"),
        }
    }
}

#[test]
fn original_coulomb_law_is_retained_and_reevaluated_after_heating() {
    let state = frozen_state(0.001, 6000.0, 12000.0, 1.0);
    let (base, card) = source_with_card(&state);
    let original = base
        .with_coulomb_drag(coulomb_model(&state), || false)
        .unwrap();
    let next = heated(&state, 6500.0, 15000.0);
    let refreshed = original.at_frozen_state(&next, &[&card], || false).unwrap();
    assert_eq!(refreshed.coulomb_model(), original.coulomb_model());
    assert_ne!(refreshed.coulomb_pairs(), original.coulomb_pairs());
    let independently_rebuilt = source_with_card(&next)
        .0
        .with_coulomb_drag(original.coulomb_model().unwrap().clone(), || false)
        .unwrap();
    relative(
        refreshed.conductivity_s_m(),
        independently_rebuilt.conductivity_s_m(),
        2e-12,
    );
}

#[test]
fn refresh_refuses_changed_populations_volume_cards_source_domains_and_cancellation() {
    let state = frozen_state(0.001, 6000.0, 12000.0, 1.0);
    let (original, card) = source_with_card(&state);
    let saved = original.clone();
    let next = heated(&state, 6500.0, 15000.0);
    assert!(original.at_frozen_state(&next, &[], || false).is_err());
    assert!(
        original
            .at_frozen_state(&next, &[&card, &card], || false)
            .is_err()
    );
    let (_, wrong) = collision_source_and_card(
        "H+",
        &[("T_h", 6000.0), ("T_e", 12000.0)],
        state.pressure().value(),
        2e9,
    );
    assert!(
        original
            .at_frozen_state(&next, &[&wrong], || false)
            .is_err()
    );
    assert!(
        original
            .at_frozen_state(&next, &[&card, &wrong], || false)
            .is_err()
    );
    let changed = frozen_state(0.002, 6500.0, 15000.0, 1.0);
    assert!(
        original
            .at_frozen_state(&changed, &[&card], || false)
            .is_err()
    );
    let changed_volume = FrozenIonizedGasStateV1::new(
        state.thermodynamics(),
        state.amounts(),
        Volume::new(2.0),
        Temperature::new(6500.0),
        Temperature::new(15000.0),
        1e-12,
        || false,
    )
    .unwrap();
    assert!(
        original
            .at_frozen_state(&changed_volume, &[&card], || false)
            .is_err()
    );
    assert_eq!(
        original.at_frozen_state(&next, &[&card], || true),
        Err(IonizedGasRefreshError::Cancelled)
    );
    let mut polls = 0;
    assert_eq!(
        original.at_frozen_state(&next, &[&card], || {
            polls += 1;
            polls > 3
        }),
        Err(IonizedGasRefreshError::Cancelled)
    );
    assert_eq!(original, saved);
}

#[test]
fn common_temperature_source_remains_common_temperature_only() {
    let state = frozen_state(0.001, 8000.0, 8000.0, 1.0);
    let (ion, ion_card) =
        collision_source_and_card("H+", &[("T", 8000.0)], state.pressure().value(), 1e9);
    let (electron, electron_card) =
        collision_source_and_card("e-", &[("T", 8000.0)], state.pressure().value(), 1e10);
    let original =
        IonizedGasConductivity::from_state(state.clone(), vec![ion, electron], limits(), || false)
            .unwrap();
    assert!(
        original
            .at_frozen_state(
                &heated(&state, 8500.0, 8500.0),
                &[&ion_card, &electron_card],
                || false
            )
            .is_ok()
    );
    assert!(
        original
            .at_frozen_state(
                &heated(&state, 8500.0, 9000.0),
                &[&ion_card, &electron_card],
                || false
            )
            .is_err()
    );
}

/// Every f64 a refresh must reproduce: the source's own values, each carrier,
/// its relaxation, and the resolved values or contributions behind that rate.
fn physical_bits(source: &IonizedGasConductivity) -> Vec<u64> {
    let mut values = vec![
        source.conductivity_s_m(),
        source.charged_mole_fraction(),
        source.coulomb_coupling(),
        source.maximum_thermal_phase_space_density(),
        source.maximum_electric_field_v_m(),
    ];
    for carrier in source.carriers() {
        let relaxation = carrier.relaxation();
        values.extend([
            carrier.number_density_per_m3(),
            carrier.mass_kg(),
            carrier.temperature_k(),
            carrier.drift_mobility_m2_v_s(),
            carrier.conductivity_s_m(),
            carrier.maximum_electric_field_v_m(),
            relaxation.rate_per_s(),
            relaxation.heavy_temperature_k(),
            relaxation.electron_temperature_k(),
            relaxation.pressure_pa(),
        ]);
        match relaxation.source() {
            CarrierMomentumSource::Material { material, .. } => {
                values.extend(material.query_point().iter().map(|(_, value)| *value));
                values.extend(material.properties().iter().map(|p| p.value_si()));
            }
            CarrierMomentumSource::ElectronNeutral {
                contributions,
                maximum_target_thermal_speed_ratio,
                ..
            } => {
                values.push(*maximum_target_thermal_speed_ratio);
                for c in contributions {
                    values.extend([
                        c.number_density_per_m3(),
                        c.rate_per_s(),
                        c.target_thermal_speed_ratio(),
                    ]);
                }
            }
        }
    }
    values.into_iter().map(f64::to_bits).collect()
}

#[test]
fn reacting_refresh_uses_changed_carrier_and_neutral_densities_without_lte_reset() {
    let state = frozen_state(0.001, 6000.0, 12000.0, 1.0);
    let (original, card) = source_with_card(&state);
    // Independent manufactured endpoint: population and temperature changes
    // have distinct analytic effects on constant-cross-section electron drag.
    let next = frozen_state(0.003, 6500.0, 15000.0, 1.0);
    let refreshed = original
        .at_reacting_state(&next, &[&card], || false)
        .unwrap();
    assert_eq!(refreshed.state(), &next);
    assert!(refreshed.equilibrium().is_none());
    assert_eq!(refreshed.limits(), original.limits());
    let before = original
        .carriers()
        .iter()
        .find(|s| s.is_electron())
        .unwrap();
    let after = refreshed
        .carriers()
        .iter()
        .find(|s| s.is_electron())
        .unwrap();
    let drag_ratio = (0.997 / 0.999) * (15000.0_f64 / 12000.0).sqrt();
    relative(
        after.number_density_per_m3() / before.number_density_per_m3(),
        3.0,
        2e-15,
    );
    relative(
        after.relaxation().rate_per_s() / before.relaxation().rate_per_s(),
        drag_ratio,
        2e-8,
    );
    relative(
        after.conductivity_s_m() / before.conductivity_s_m(),
        3.0 / drag_ratio,
        2e-8,
    );
    match (before.relaxation().source(), after.relaxation().source()) {
        (
            CarrierMomentumSource::ElectronNeutral {
                contributions: old, ..
            },
            CarrierMomentumSource::ElectronNeutral {
                state: bound,
                contributions: new,
                ..
            },
        ) => {
            assert_eq!(bound, &next);
            assert_eq!(old[0].evaluation().source(), new[0].evaluation().source());
            assert_eq!(
                old[0].evaluation().controls(),
                new[0].evaluation().controls()
            );
            relative(
                new[0].number_density_per_m3() / old[0].number_density_per_m3(),
                0.997 / 0.999,
                2e-15,
            );
        }
        _ => panic!("reaction refresh changed the momentum source"),
    }
    // Recombination can reduce the populations with the identical original
    // source. Returning to the old physical state recovers its transport.
    let recovered = refreshed
        .at_reacting_state(&state, &[&card], || false)
        .unwrap();
    // A refresh pins each selected claim; `original` was resolved with
    // SingleClaimOnly. The library states the rule in
    // src/liquid/density/tests.rs:283-287: "A refresh records pinned-claim
    // selection; an independent resolution records SingleClaimOnly. Their
    // physical tuples, not policy receipts, must agree exactly." So the whole
    // struct equals the original refreshed through the same pinned path, and
    // every physical value equals the single-claim original bit for bit; only
    // the receipt policy, its provenance and the identity binding them differ.
    assert_eq!(
        recovered,
        original
            .at_reacting_state(&state, &[&card], || false)
            .unwrap()
    );
    assert_eq!(recovered.state(), original.state());
    assert_eq!(physical_bits(&recovered), physical_bits(&original));
    assert!(original.at_frozen_state(&next, &[&card], || false).is_err());
}

#[test]
fn reacting_refresh_rebinds_original_coulomb_law_to_actual_populations() {
    let state = frozen_state(0.001, 6000.0, 12000.0, 1.0);
    let (base, card) = source_with_card(&state);
    let original = base
        .with_coulomb_drag(coulomb_model(&state), || false)
        .unwrap();
    let next = frozen_state(0.002, 6500.0, 15000.0, 1.0);
    let refreshed = original
        .at_reacting_state(&next, &[&card], || false)
        .unwrap();
    let old_law = original.coulomb_model().unwrap();
    let new_law = refreshed.coulomb_model().unwrap();
    assert!(new_law.is_bound_to_populations(&next));
    assert!(!old_law.is_bound_to_populations(&next));
    assert_eq!(new_law.screening(), old_law.screening());
    assert_eq!(new_law.cutoff(), old_law.cutoff());
    assert_eq!(new_law.limits(), old_law.limits());
    assert_ne!(refreshed.coulomb_pairs(), original.coulomb_pairs());
    let independent = source_with_card(&next)
        .0
        .with_coulomb_drag(coulomb_model(&next), || false)
        .unwrap();
    relative(
        refreshed.conductivity_s_m(),
        independent.conductivity_s_m(),
        2e-12,
    );
    let field = 0.01 * refreshed.maximum_electric_field_v_m();
    relative(
        refreshed
            .joule_heating_w_m3([field, 0.0])
            .unwrap()
            .total_power_density_w_m3(),
        independent
            .joule_heating_w_m3([field, 0.0])
            .unwrap()
            .total_power_density_w_m3(),
        2e-12,
    );
}

#[test]
fn reacting_refresh_preserves_source_support_domains_and_atomic_publication() {
    let state = frozen_state(0.001, 6000.0, 12000.0, 1.0);
    let (original, card) = source_with_card(&state);
    let saved = original.clone();
    let next = frozen_state(0.002, 6500.0, 15000.0, 1.0);
    assert!(original.at_reacting_state(&next, &[], || false).is_err());
    assert!(
        original
            .at_reacting_state(&next, &[&card, &card], || false)
            .is_err()
    );
    for unsupported in [
        frozen_state(0.0, 6500.0, 15000.0, 1.0),
        frozen_state(0.002, 6500.0, 15000.0, 2.0),
        frozen_state(0.5, 6500.0, 15000.0, 1.0),
    ] {
        assert!(
            original
                .at_reacting_state(&unsupported, &[&card], || false)
                .is_err()
        );
    }
    assert_eq!(
        original.at_reacting_state(&next, &[&card], || true),
        Err(IonizedGasRefreshError::Cancelled)
    );
    let mut polls = 0;
    assert_eq!(
        original.at_reacting_state(&next, &[&card], || {
            polls += 1;
            polls > 8
        }),
        Err(IonizedGasRefreshError::Cancelled)
    );
    assert_eq!(original, saved);
}
