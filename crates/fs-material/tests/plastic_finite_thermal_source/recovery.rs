//! Original-source identity and complete recovery state on finite thermal laws.

use super::*;
use fs_material::plastic::finite::thermal::FiniteThermalInelasticMaterial;
use fs_material::plastic::finite::thermal::source::recovery::*;

const SELECTORS: [(&str, &str); 3] = [
    (
        THERMAL_FINITE_RECOVERING_NORTON_LAW_ID,
        THERMAL_FINITE_PRIMARY_NORTON_LAW_ID,
    ),
    (
        THERMAL_FINITE_RECOVERING_NORTON_J2_LAW_ID,
        THERMAL_FINITE_PRIMARY_NORTON_J2_LAW_ID,
    ),
    (
        THERMAL_FINITE_RECOVERING_NORTON_J2_VOCE_LAW_ID,
        THERMAL_FINITE_PRIMARY_NORTON_J2_VOCE_LAW_ID,
    ),
];

fn recovering_card(selector: &str, base_selector: &str, rate: f64) -> ConstitutiveModelCard {
    let mut card = primary_fixture(base_selector);
    card.law = LawId(selector.into());
    card.state_schema_version = 3;
    for (name, value, dims) in [
        (
            "creep_reference_recovery_rate_s",
            rate,
            Dims([0, 0, -1, 0, 0, 0]),
        ),
        ("creep_recovery_activation_temperature_k", 1800.0, K),
    ] {
        card.parameters
            .insert(name.into(), LawParameter { value, dims });
    }
    card
}

#[test]
fn original_source_and_complete_recovery_history_survive_finite_heated_flow() {
    let identity = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
    let loaded = [1.08, 0.0, 0.0, 0.0, 0.97, 0.0, 0.0, 0.0, 0.99];
    for (selector, base_selector) in SELECTORS {
        let card = recovering_card(selector, base_selector, 0.4);
        let source = ResolvedAffineThermalRecoveringCreep::try_from_model_card_at(
            &card,
            &context(),
            "temperature",
        )
        .unwrap();
        assert_eq!(source.model_card(), &card);
        assert_eq!(source.context(), &context());
        assert_eq!(source.temperature_axis(), "temperature");
        assert_eq!(source.law().recovery().activation_temperature_k, 1800.0);
        let old = source.initial_state();
        let step = source
            .checked_increment_for_step(&identity, &old, T_REF, &loaded, 345.0, 1e-5, 0.8)
            .unwrap();
        let history = &step.mechanical.state;
        let creep = if selector == THERMAL_FINITE_RECOVERING_NORTON_LAW_ID {
            history.native.alpha
        } else {
            assert!(history.native.alpha > 0.0);
            history.native.creep_alpha
        };
        let rate = 0.4 * (1800.0_f64 * (1.0 / T_REF - 1.0 / 345.0)).exp();
        assert!(creep > 0.0 && history.resistance_strain > 0.0);
        assert!((history.resistance_strain - creep / (1.0 + rate * 0.8)).abs() < 2e-14);
        assert!(step.mechanical.inelastic_dissipation_j_m3() > 0.0);
        assert!(step.balance_error_j_m3.abs() < 1e-7);
        assert_eq!(old, source.initial_state());

        // Explicit zero recovery preserves the old primary solution exactly,
        // while retaining the new source identity and full history schema.
        let zero_card = recovering_card(selector, base_selector, 0.0);
        let zero = ResolvedAffineThermalRecoveringCreep::try_from_model_card_at(
            &zero_card,
            &context(),
            "temperature",
        )
        .unwrap();
        let reference = ResolvedAffineThermalFiniteJ2::try_from_model_card_at(
            &primary_fixture(base_selector),
            &context(),
            "temperature",
        )
        .unwrap();
        let expected = reference
            .law()
            .checked_increment_for_step(
                &identity,
                &reference.law().initial_state(),
                T_REF,
                &loaded,
                345.0,
                1e-5,
                0.8,
            )
            .unwrap();
        let actual = zero
            .checked_increment_for_step(
                &identity,
                &zero.initial_state(),
                T_REF,
                &loaded,
                345.0,
                1e-5,
                0.8,
            )
            .unwrap();
        assert_eq!(actual.mechanical.state.native, expected.mechanical.state);
        assert_eq!(
            actual.mechanical.first_piola_pa,
            expected.mechanical.first_piola_pa
        );
        assert_eq!(
            actual.mechanical.kirchhoff_pa,
            expected.mechanical.kirchhoff_pa
        );
        assert_eq!(actual.mechanical.tangent_pa, expected.mechanical.tangent_pa);
        assert_eq!(
            actual.mechanical.incremental_potential_j_m3,
            expected.mechanical.incremental_potential_j_m3
        );
        assert_eq!(
            actual.mechanical.flow_potential_j_m3,
            expected.mechanical.flow_potential_j_m3
        );
        assert_eq!(
            actual.mechanical.stored_energy_j_m3,
            expected.mechanical.stored_energy_j_m3
        );
        assert_eq!(
            actual.mechanical.plastic_dissipation_j_m3,
            expected.mechanical.plastic_dissipation_j_m3
        );
        assert_eq!(
            actual.mechanical.creep_dissipation_j_m3,
            expected.mechanical.creep_dissipation_j_m3
        );
        assert_eq!(actual.reversible_heat_j_m3, expected.reversible_heat_j_m3);
        assert_eq!(
            actual.thermoelastic_storage_change_j_m3,
            expected.thermoelastic_storage_change_j_m3
        );
    }
}

#[test]
fn recovery_requires_exact_law_state_units_and_temperature_domain() {
    let (selector, base_selector) = SELECTORS[1];
    let card = recovering_card(selector, base_selector, 0.4);
    let rejects = |candidate: &ConstitutiveModelCard| {
        assert!(
            ResolvedAffineThermalRecoveringCreep::try_from_model_card_at(
                candidate,
                &context(),
                "temperature",
            )
            .is_err()
        );
    };
    let mut changed = card.clone();
    changed.state_schema_version = 2;
    rejects(&changed);
    changed = card.clone();
    changed
        .parameters
        .get_mut("creep_reference_recovery_rate_s")
        .unwrap()
        .dims = Dims::NONE;
    rejects(&changed);
    changed = card.clone();
    changed
        .parameters
        .get_mut("creep_recovery_activation_temperature_k")
        .unwrap()
        .value = -1.0;
    rejects(&changed);
    changed = card.clone();
    changed
        .parameters
        .remove("creep_hardening_reference_strain");
    rejects(&changed);
    changed = card.clone();
    changed.validity = ValidityDomain::unconstrained().with_quantity(
        "temperature",
        semantic(QuantityKind::AbsoluteTemperature),
        300.0,
        350.0,
    );
    rejects(&changed);
    rejects(&primary_fixture(base_selector));
    assert!(
        ResolvedAffineThermalFiniteJ2::try_from_model_card_at(&card, &context(), "temperature",)
            .is_err()
    );
}
