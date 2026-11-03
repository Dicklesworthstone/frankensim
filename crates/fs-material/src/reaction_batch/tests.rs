use super::*;
use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{LawId, LawParameter, Provenance};
use fs_qty::Frequency;
use fs_thermochem::stoichiometry::minerals::MineralReaction;
use std::collections::BTreeMap;

fn source(k: f64) -> ConstitutiveModelCard {
    ConstitutiveModelCard {
        law: LawId(BILINEAR_BATCH_LAW.into()),
        law_version: 1,
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        parameters: BTreeMap::from([(
            RATE_PARAMETER.into(),
            LawParameter {
                value: k,
                dims: BilinearRateConstant::DIMS,
            },
        )]),
        validity: ValidityDomain::unconstrained()
            .with_quantity(
                ELAPSED_TIME_AXIS,
                QuantitySpec::dimensional(Time::DIMS),
                0.0,
                1e30,
            )
            .with_quantity(
                VOLUME_AXIS,
                QuantitySpec::dimensional(Volume::DIMS),
                1e-6,
                1e6,
            )
            .with_quantity(
                TEMPERATURE_AXIS,
                QuantitySpec::dimensional(Temperature::DIMS),
                273.15,
                373.15,
            )
            .with_quantity(
                REACTANT_A_CONCENTRATION_AXIS,
                QuantitySpec::dimensional(AmountConcentration::DIMS),
                0.0,
                1e12,
            )
            .with_quantity(
                REACTANT_B_CONCENTRATION_AXIS,
                QuantitySpec::dimensional(AmountConcentration::DIMS),
                0.0,
                1e12,
            ),
        sources: vec![hash_domain(
            "bilinear-batch-test",
            b"authored mathematical fixture, not a measured mineral rate",
        )],
        provenance: Provenance {
            source: "authored closed-form rate-law fixture; not calibrated".into(),
            license: "project fixture".into(),
            artifact: None,
        },
    }
}

fn model(kind: MineralReaction, mineral: &str, k: f64) -> BilinearBatchReaction {
    BilinearBatchReaction::from_card(
        kind.build().unwrap(),
        SpeciesId::new(mineral).unwrap(),
        SpeciesId::new("H+(aq)").unwrap(),
        source(k),
    )
    .unwrap()
}

fn calcite(k: f64) -> BilinearBatchReaction {
    model(
        MineralReaction::CalciteBicarbonateNeutralization,
        "CaCO3(s)",
        k,
    )
}

fn inventory(model: &BilinearBatchReaction, a: f64, b: f64) -> Vec<(SpeciesId, Amount)> {
    model
        .reaction()
        .terms()
        .iter()
        .map(|term| {
            let amount = if term.species == model.reactants()[0] {
                a
            } else if term.species == model.reactants()[1] {
                b
            } else {
                0.0
            };
            (term.species.clone(), Amount::new(amount))
        })
        .collect()
}

fn advance(
    model: &BilinearBatchReaction,
    initial: &[(SpeciesId, Amount)],
    t: f64,
) -> BatchReactionAdvance {
    model
        .advance(
            initial,
            Volume::new(1.0),
            Time::new(t),
            Temperature::new(298.15),
            &QueryPoint::new(),
        )
        .unwrap()
}

fn amount(outcome: &BatchReactionAdvance, species: &str) -> f64 {
    stock(&outcome.final_inventory, &SpeciesId::new(species).unwrap())
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

#[test]
fn equal_capacity_known_answer_retains_equation_and_original_source() {
    let model = calcite(0.01);
    let mut initial = inventory(&model, 10.0, 10.0);
    initial.reverse(); // Input order does not become species identity.
    let result = advance(&model, &initial, 10.0);
    // dξ/dt=.01(10-ξ)^2: at t=10, 1/(10-ξ)=1/10+.01*10.
    close(result.extent.value(), 5.0, 2e-14);
    close(amount(&result, "CaCO3(s)"), 5.0, 2e-14);
    close(amount(&result, "H+(aq)"), 5.0, 2e-14);
    close(amount(&result, "Ca[2+](aq)"), 5.0, 2e-14);
    close(amount(&result, "HCO3[-](aq)"), 5.0, 2e-14);
    assert_eq!(result.limiting_extent.value(), 10.0);
    assert_eq!(result.model().source(), model.source());
    assert_eq!(result.model().reaction().id(), model.reaction().id());
    assert_eq!(result.model().reaction().terms(), model.reaction().terms());
    assert_eq!(result.model().reactants(), model.reactants());
    assert_eq!(result.initial_point.axes()[ELAPSED_TIME_AXIS], 0.0);
    assert_eq!(result.final_point.axes()[ELAPSED_TIME_AXIS], 10.0);
    assert_eq!(
        result.initial_point.axes()[REACTANT_A_CONCENTRATION_AXIS],
        10.0
    );
    close(
        result.final_point.axes()[REACTANT_A_CONCENTRATION_AXIS],
        5.0,
        2e-14,
    );
    assert_eq!(result.initial_inventory, inventory(&model, 10.0, 10.0));
    assert_eq!(result.amounts.len(), model.reaction().terms().len());
    for account in &result.reactants {
        close(
            account.initial.value(),
            account.consumed.value() + account.remaining.value(),
            2e-14,
        );
    }
}

#[test]
fn unequal_capacity_known_answer_and_role_reversal() {
    let model = calcite(0.01);
    let t = core::f64::consts::LN_2 / 0.1;
    // Integration gives (10-ξ)/(20-ξ)=1/4, hence ξ=20/3.
    for (a, b) in [(10.0, 20.0), (20.0, 10.0)] {
        let result = advance(&model, &inventory(&model, a, b), t);
        close(result.extent.value(), 20.0 / 3.0, 3e-14);
        close(result.reactants[0].remaining.value(), a - 20.0 / 3.0, 3e-14);
        close(result.reactants[1].remaining.value(), b - 20.0 / 3.0, 3e-14);
    }
    let nearly_equal = advance(
        &model,
        &inventory(&model, 10.0, 10.0 * (1.0 + f64::EPSILON)),
        10.0,
    );
    close(nearly_equal.extent.value(), 5.0, 5e-14);
}

#[test]
fn dolomite_uses_two_acid_moles_and_depletes_finite_reagent_during_progress() {
    let model = model(
        MineralReaction::DolomiteBicarbonateNeutralization,
        "CaMg(CO3)2(s)",
        0.005,
    );
    // dξ/dt=.005(10-ξ)(20-2ξ)=.01(10-ξ)^2.
    let result = advance(&model, &inventory(&model, 10.0, 20.0), 10.0);
    close(result.extent.value(), 5.0, 2e-14);
    close(result.reactants[0].consumed.value(), 5.0, 2e-14);
    close(result.reactants[1].consumed.value(), 10.0, 4e-14);
    close(result.reactants[1].remaining.value(), 10.0, 4e-14);
    assert_eq!(result.reactants[1].coefficient, 2);
    close(amount(&result, "Ca[2+](aq)"), 5.0, 2e-14);
    close(amount(&result, "Mg[2+](aq)"), 5.0, 2e-14);
    close(amount(&result, "HCO3[-](aq)"), 10.0, 4e-14);
    let acid_limited = advance(&model, &inventory(&model, 10.0, 2.0), 100.0);
    assert_eq!(acid_limited.limiting_extent.value(), 1.0);
    assert!(acid_limited.extent.value() > 0.99 && acid_limited.extent.value() < 1.0);
    assert!(acid_limited.reactants[1].remaining.value() > 0.0);
    assert!(acid_limited.reactants[0].remaining.value() > 9.0);
    // With the same mineral, time and volume, lower acid also lowers early rate.
    let low = advance(&model, &inventory(&model, 10.0, 2.0), 1.0);
    let high = advance(&model, &inventory(&model, 10.0, 20.0), 1.0);
    assert!(low.extent.value() < high.extent.value() / 5.0);
}

#[test]
fn sequential_contact_matches_one_batch_and_preserves_initial_products() {
    let model = model(
        MineralReaction::DolomiteBicarbonateNeutralization,
        "CaMg(CO3)2(s)",
        0.005,
    );
    let mut initial = inventory(&model, 12.0, 20.0);
    initial
        .iter_mut()
        .find(|(id, _)| id.as_str() == "Mg[2+](aq)")
        .unwrap()
        .1 = Amount::new(3.0);
    let one = advance(&model, &initial, 10.0);
    let first = advance(&model, &initial, 6.0);
    let second = advance(&model, &first.final_inventory, 4.0);
    close(
        first.extent.value() + second.extent.value(),
        one.extent.value(),
        8e-14,
    );
    for ((a, actual), (b, expected)) in second.final_inventory.iter().zip(&one.final_inventory) {
        assert_eq!(a, b);
        close(actual.value(), expected.value(), 2e-13);
    }
    close(amount(&one, "Mg[2+](aq)"), 3.0 + one.extent.value(), 2e-14);
}

#[test]
fn actual_volume_controls_rate_and_both_concentration_coordinates() {
    let model = calcite(0.01);
    let initial = inventory(&model, 10.0, 10.0);
    let result = model
        .advance(
            &initial,
            Volume::new(2.0),
            Time::new(10.0),
            Temperature::new(298.15),
            &QueryPoint::new(),
        )
        .unwrap();
    // At V=2, 1/(10-ξ)=1/10+.005*10, so ξ=10/3.
    close(result.extent.value(), 10.0 / 3.0, 3e-14);
    assert_eq!(
        result.initial_point.axes()[REACTANT_A_CONCENTRATION_AXIS],
        5.0
    );
    assert_eq!(
        result.initial_point.axes()[REACTANT_B_CONCENTRATION_AXIS],
        5.0
    );
    close(
        result.final_point.axes()[REACTANT_A_CONCENTRATION_AXIS],
        10.0 / 3.0,
        3e-14,
    );
}

#[test]
fn zero_time_rate_or_reactant_returns_complete_unchanged_inventory() {
    for (k, a, b, time) in [
        (0.01, 10.0, 10.0, 0.0),
        (0.0, 10.0, 10.0, 10.0),
        (0.01, 0.0, 10.0, 10.0),
        (0.01, 10.0, 0.0, 10.0),
        (0.01, 0.0, 0.0, 10.0),
    ] {
        let model = calcite(k);
        let initial = inventory(&model, a, b);
        let result = advance(&model, &initial, time);
        assert_eq!(result.extent.value(), 0.0);
        assert_eq!(result.final_inventory, initial);
        assert!(result.amounts.iter().all(|flow| flow.amount.value() == 0.0));
    }
}

#[test]
fn source_must_cover_actual_final_concentrations_and_all_fixed_conditions() {
    let base = calcite(0.01);
    let initial = inventory(&base, 10.0, 10.0);
    let mut card = source(0.01);
    card.validity = card.validity.with_quantity(
        REACTANT_A_CONCENTRATION_AXIS,
        QuantitySpec::dimensional(AmountConcentration::DIMS),
        8.0,
        12.0,
    );
    let restricted = BilinearBatchReaction::from_card(
        base.reaction().clone(),
        base.reactants()[0].clone(),
        base.reactants()[1].clone(),
        card,
    )
    .unwrap();
    assert_eq!(
        restricted
            .advance(
                &initial,
                Volume::new(1.0),
                Time::new(10.0),
                Temperature::new(298.15),
                &QueryPoint::new()
            )
            .unwrap_err(),
        BatchReactionError::OutsideValidity { endpoint: "final" }
    );

    let mut card = source(0.01);
    card.validity = card.validity.with_quantity(
        "agitation-frequency",
        QuantitySpec::dimensional(Frequency::DIMS),
        1.0,
        100.0,
    );
    let restricted = BilinearBatchReaction::from_card(
        base.reaction().clone(),
        base.reactants()[0].clone(),
        base.reactants()[1].clone(),
        card,
    )
    .unwrap();
    // Even zero-time calculations cannot bypass source admission.
    assert_eq!(
        restricted
            .advance(
                &initial,
                Volume::new(1.0),
                Time::new(0.0),
                Temperature::new(298.15),
                &QueryPoint::new()
            )
            .unwrap_err(),
        BatchReactionError::OutsideValidity {
            endpoint: "initial"
        }
    );
    let mistyped = QueryPoint::new()
        .with_quantity(
            "agitation-frequency",
            QuantitySpec::dimensional(Time::DIMS),
            30.0,
        )
        .unwrap();
    assert!(matches!(
        restricted.advance(
            &initial,
            Volume::new(1.0),
            Time::new(10.0),
            Temperature::new(298.15),
            &mistyped
        ),
        Err(BatchReactionError::OutsideValidity { .. })
    ));
    let correct = QueryPoint::new()
        .with_quantity(
            "agitation-frequency",
            QuantitySpec::dimensional(Frequency::DIMS),
            30.0,
        )
        .unwrap();
    assert!(
        restricted
            .advance(
                &initial,
                Volume::new(1.0),
                Time::new(10.0),
                Temperature::new(298.15),
                &correct
            )
            .is_ok()
    );
}

#[test]
fn owned_axes_cannot_be_spoofed_and_invalid_time_temperature_volume_refuse() {
    let model = calcite(0.01);
    let initial = inventory(&model, 10.0, 10.0);
    for axis in OWNED_AXES {
        let spoofed = QueryPoint::new().with(axis, 1.0).unwrap();
        assert!(matches!(
            model.advance(
                &initial,
                Volume::new(1.0),
                Time::new(10.0),
                Temperature::new(298.15),
                &spoofed
            ),
            Err(BatchReactionError::Invalid(_))
        ));
    }
    for (v, t, temp) in [
        (0.0, 10.0, 298.15),
        (-1.0, 10.0, 298.15),
        (1.0, -1.0, 298.15),
        (1.0, 10.0, 0.0),
        (1.0, f64::INFINITY, 298.15),
        (1.0, 10.0, f64::NAN),
    ] {
        assert!(matches!(
            model.advance(
                &initial,
                Volume::new(v),
                Time::new(t),
                Temperature::new(temp),
                &QueryPoint::new()
            ),
            Err(BatchReactionError::Invalid(_))
        ));
    }
    assert!(matches!(
        model.advance(
            &initial,
            Volume::new(1.0),
            Time::new(10.0),
            Temperature::new(500.0),
            &QueryPoint::new()
        ),
        Err(BatchReactionError::OutsideValidity { .. })
    ));
    assert_eq!(
        model
            .advance(
                &initial,
                Volume::new(1.0),
                Time::new(2e30),
                Temperature::new(298.15),
                &QueryPoint::new()
            )
            .unwrap_err(),
        BatchReactionError::OutsideValidity { endpoint: "final" }
    );
}

#[test]
fn unsupported_laws_sources_and_reactant_mappings_refuse() {
    let base = calcite(0.01);
    for case in 0..6 {
        let mut card = source(0.01);
        match case {
            0 => card.sources.clear(),
            1 => card.parameters.get_mut(RATE_PARAMETER).unwrap().dims = Frequency::DIMS,
            2 => card.parameters.get_mut(RATE_PARAMETER).unwrap().value = -0.1,
            3 => card.provenance.source.clear(),
            4 => {
                card.validity = card.validity.with_quantity(
                    ELAPSED_TIME_AXIS,
                    QuantitySpec::dimensional(Time::DIMS),
                    1.0,
                    100.0,
                )
            }
            _ => {
                card.validity = card.validity.with_quantity(
                    REACTANT_A_CONCENTRATION_AXIS,
                    QuantitySpec::dimensional(Amount::DIMS),
                    0.0,
                    10.0,
                )
            }
        }
        assert!(
            BilinearBatchReaction::from_card(
                base.reaction().clone(),
                base.reactants()[0].clone(),
                base.reactants()[1].clone(),
                card
            )
            .is_err()
        );
    }
    assert!(
        BilinearBatchReaction::from_card(
            base.reaction().clone(),
            base.reactants()[0].clone(),
            base.reactants()[0].clone(),
            source(0.01)
        )
        .is_err()
    );
    assert!(
        BilinearBatchReaction::from_card(
            base.reaction().clone(),
            SpeciesId::new("Ca[2+](aq)").unwrap(),
            base.reactants()[1].clone(),
            source(0.01)
        )
        .is_err()
    );
    assert!(
        BilinearBatchReaction::from_card(
            MineralReaction::GoldCyanidation.build().unwrap(),
            base.reactants()[0].clone(),
            base.reactants()[1].clone(),
            source(0.01)
        )
        .is_err()
    );
}

#[test]
fn incomplete_inventory_and_unrepresentable_progress_refuse_without_mutation() {
    let model = calcite(0.01);
    let initial = inventory(&model, 10.0, 10.0);
    let saved = initial.clone();
    assert!(matches!(
        model.advance(
            &initial[..initial.len() - 1],
            Volume::new(1.0),
            Time::new(10.0),
            Temperature::new(298.15),
            &QueryPoint::new()
        ),
        Err(BatchReactionError::Reaction(_))
    ));
    // Both are finite and inside support: tiny progress cannot change large
    // stock; huge progress cannot round a still-positive remainder to zero.
    for t in [1e-30, 1e20] {
        assert!(matches!(
            model.advance(
                &initial,
                Volume::new(1.0),
                Time::new(t),
                Temperature::new(298.15),
                &QueryPoint::new()
            ),
            Err(BatchReactionError::Arithmetic(_))
                | Err(BatchReactionError::Reaction(
                    ReactionBalanceError::Arithmetic(_)
                ))
        ));
    }
    assert_eq!(initial, saved);
}
