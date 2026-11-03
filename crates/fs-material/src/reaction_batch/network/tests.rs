use super::*;
use crate::reaction_batch::{
    BILINEAR_BATCH_LAW, BilinearRateConstant, ELAPSED_TIME_AXIS, RATE_PARAMETER,
    REACTANT_A_CONCENTRATION_AXIS, REACTANT_B_CONCENTRATION_AXIS, TEMPERATURE_AXIS, VOLUME_AXIS,
};
use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_qty::{AmountConcentration, QuantitySpec};
use fs_thermochem::stoichiometry::minerals::MineralReaction;

fn card(k: f64) -> ConstitutiveModelCard {
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
                100.0,
            )
            .with_quantity(
                VOLUME_AXIS,
                QuantitySpec::dimensional(Volume::DIMS),
                1.0,
                1.0,
            )
            .with_quantity(
                TEMPERATURE_AXIS,
                QuantitySpec::dimensional(Temperature::DIMS),
                298.15,
                298.15,
            )
            .with_quantity(
                REACTANT_A_CONCENTRATION_AXIS,
                QuantitySpec::dimensional(AmountConcentration::DIMS),
                0.0,
                100.0,
            )
            .with_quantity(
                REACTANT_B_CONCENTRATION_AXIS,
                QuantitySpec::dimensional(AmountConcentration::DIMS),
                0.0,
                100.0,
            ),
        sources: vec![hash_domain(
            "simultaneous-batch-fixture",
            b"authored bilinear constants; not measured mineral kinetics",
        )],
        provenance: Provenance {
            source: "synthetic competition fixture".into(),
            license: "project fixture".into(),
            artifact: None,
        },
    }
}
fn copper(k: f64) -> BilinearBatchReaction {
    BilinearBatchReaction::from_card(
        MineralReaction::CupricOxideAcidLeach.build().unwrap(),
        SpeciesId::new("CuO(s)").unwrap(),
        SpeciesId::new("H+(aq)").unwrap(),
        card(k),
    )
    .unwrap()
}
fn dolomite(k: f64) -> BilinearBatchReaction {
    BilinearBatchReaction::from_card(
        MineralReaction::DolomiteBicarbonateNeutralization
            .build()
            .unwrap(),
        SpeciesId::new("CaMg(CO3)2(s)").unwrap(),
        SpeciesId::new("H+(aq)").unwrap(),
        card(k),
    )
    .unwrap()
}
fn network(models: Vec<BilinearBatchReaction>) -> BilinearBatchNetwork {
    BilinearBatchNetwork::new(models, SpeciesId::new("H+(aq)").unwrap()).unwrap()
}
fn inventory(network: &BilinearBatchNetwork, acid: f64) -> Vec<(SpeciesId, Amount)> {
    network
        .reaction_network()
        .species()
        .iter()
        .map(|id| {
            let value = match id.as_str() {
                "CuO(s)" | "CaMg(CO3)2(s)" => 10.0,
                "H+(aq)" => acid,
                _ => 0.0,
            };
            (id.clone(), Amount::new(value))
        })
        .collect()
}
fn conditions(network: &BilinearBatchNetwork) -> Vec<(ReactionId, QueryPoint)> {
    network
        .models()
        .iter()
        .map(|model| (model.reaction().id().clone(), QueryPoint::new()))
        .collect()
}
fn budget() -> BatchIntegrationBudget {
    BatchIntegrationBudget {
        maximum_step: Time::new(1.0),
        absolute_extent_error: Amount::new(1e-9),
        relative_extent_error: Dimensionless::new(1e-8),
        maximum_attempts: 4096,
    }
}
fn run(network: &BilinearBatchNetwork, acid: f64, seconds: f64) -> BatchNetworkAdvance {
    network
        .advance(
            &inventory(network, acid),
            Volume::new(1.0),
            Time::new(seconds),
            Temperature::new(298.15),
            &conditions(network),
            budget(),
        )
        .unwrap()
}
fn amount(result: &BatchNetworkAdvance, species: &str) -> f64 {
    result
        .stage
        .final_inventory
        .iter()
        .find(|(id, _)| id.as_str() == species)
        .unwrap()
        .1
        .value()
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

#[test]
fn single_model_agrees_with_existing_closed_form_and_retains_original_evidence() {
    let model = copper(0.005);
    let network = network(vec![model.clone()]);
    let initial = inventory(&network, 20.0);
    let exact = model
        .advance(
            &initial,
            Volume::new(1.0),
            Time::new(10.0),
            Temperature::new(298.15),
            &QueryPoint::new(),
        )
        .unwrap();
    let result = run(&network, 20.0, 10.0);
    // n_H=2*n_CuO; dn_CuO/dt=-.01*n_CuO^2, hence n_CuO(10)=5.
    close(result.stage.reactions[0].extent.value(), 5.0, 2e-7);
    close(
        result.stage.reactions[0].extent.value(),
        exact.extent.value(),
        2e-7,
    );
    for ((id, actual), (other, expected)) in result
        .stage
        .final_inventory
        .iter()
        .zip(exact.final_inventory)
    {
        assert_eq!(id, &other);
        close(actual.value(), expected.value(), 4e-7);
    }
    assert_eq!(result.model().models()[0].source(), model.source());
    assert_eq!(
        result.model().models()[0].reaction().terms(),
        model.reaction().terms()
    );
    assert_eq!(
        result.source_paths[0].initial.axes()[ELAPSED_TIME_AXIS],
        0.0
    );
    assert_eq!(
        result.source_paths[0].final_point.axes()[ELAPSED_TIME_AXIS],
        10.0
    );
    assert!(result.accepted_steps > 1);
    assert!(result.estimated_extent_errors[0].1.value() <= 1.01e-7);
}

#[test]
fn symmetric_competition_matches_independent_aggregate_solution() {
    let network = network(vec![copper(0.005), dolomite(0.005)]);
    let result = run(&network, 40.0, 10.0);
    // Equal minerals stay equal: n_H=4*n_mineral. Thus dn/dt=-.02*n^2,
    // n(10)=10/3 and each reaction extent is 20/3, not an independent batch.
    for reaction in &result.stage.reactions {
        close(reaction.extent.value(), 20.0 / 3.0, 3e-7);
    }
    close(amount(&result, "CuO(s)"), 10.0 / 3.0, 3e-7);
    close(amount(&result, "CaMg(CO3)2(s)"), 10.0 / 3.0, 3e-7);
    close(amount(&result, "H+(aq)"), 40.0 / 3.0, 1.2e-6);
}

#[test]
fn reaction_inventory_and_condition_permutations_do_not_change_results() {
    let first = network(vec![copper(0.01), dolomite(0.002)]);
    let second = network(vec![dolomite(0.002), copper(0.01)]);
    let mut initial = inventory(&second, 12.0);
    let mut queries = conditions(&second);
    initial.reverse();
    queries.reverse();
    let a = run(&first, 12.0, 8.0);
    let b = second
        .advance(
            &initial,
            Volume::new(1.0),
            Time::new(8.0),
            Temperature::new(298.15),
            &queries,
            budget(),
        )
        .unwrap();
    assert_eq!(a.stage.final_inventory, b.stage.final_inventory);
    assert_eq!(a.estimated_extent_errors, b.estimated_extent_errors);
    assert_eq!(a.source_paths, b.source_paths);
    assert_eq!(a.accepted_steps, b.accepted_steps);
    assert_eq!(a.rejected_steps, b.rejected_steps);
    for (a, b) in a.stage.reactions.iter().zip(b.stage.reactions) {
        assert_eq!(a.reaction, b.reaction);
        assert_eq!(a.extent, b.extent);
    }
}

#[test]
fn unequal_proton_coefficients_share_the_same_acid_inventory() {
    let malachite = BilinearBatchReaction::from_card(
        MineralReaction::MalachiteCarbonicAcidLeach.build().unwrap(),
        SpeciesId::new("Cu2CO3(OH)2(s)").unwrap(),
        SpeciesId::new("H+(aq)").unwrap(),
        card(0.01),
    )
    .unwrap();
    let network = network(vec![copper(0.01), malachite]);
    let initial = network
        .reaction_network()
        .species()
        .iter()
        .map(|id| {
            let value = match id.as_str() {
                "CuO(s)" | "Cu2CO3(OH)2(s)" => 2.0,
                "H+(aq)" => 12.0,
                _ => 0.0,
            };
            (id.clone(), Amount::new(value))
        })
        .collect::<Vec<_>>();
    let result = network
        .advance(
            &initial,
            Volume::new(1.0),
            Time::new(10.0),
            Temperature::new(298.15),
            &conditions(&network),
            budget(),
        )
        .unwrap();
    // Equal mineral stocks and k remain equal, but consume 2+4=6 proton
    // equivalents together: H=6*n and dn/dt=-6*.01*n^2.
    let expected = 2.0 - 2.0 / (1.0 + 6.0 * 0.01 * 2.0 * 10.0);
    for reaction in &result.stage.reactions {
        close(reaction.extent.value(), expected, 1e-7);
    }
    assert_eq!(
        result.stage.reactions[0].extent,
        result.stage.reactions[1].extent
    );
    close(amount(&result, "H+(aq)"), 12.0 - 6.0 * expected, 6e-7);
    close(amount(&result, "Cu[2+](aq)"), 3.0 * expected, 3e-7);
}

#[test]
fn differing_rates_compete_for_finite_acid_without_priority_or_exhaustion() {
    let network = network(vec![copper(0.02), dolomite(0.002)]);
    let result = run(&network, 2.0, 20.0);
    let copper_remaining = amount(&result, "CuO(s)");
    let dolomite_remaining = amount(&result, "CaMg(CO3)2(s)");
    assert!(amount(&result, "H+(aq)") > 0.0);
    assert!(amount(&result, "H+(aq)") < 0.001);
    assert!(copper_remaining < dolomite_remaining);
    assert!(dolomite_remaining < 10.0); // The slower mineral still reacts.
    // Both obey d ln(n_i)/dt=-k_i*n_H/V. Their shared exposure gives this
    // independent invariant, which a sequential-priority allocation violates.
    close(
        copper_remaining / 10.0,
        (dolomite_remaining / 10.0).powi(10),
        2e-7,
    );
    let extent_sum: f64 = result
        .stage
        .reactions
        .iter()
        .map(|r| r.extent.value())
        .sum();
    assert!(extent_sum < 1.0);
    close(amount(&result, "H+(aq)"), 2.0 - 2.0 * extent_sum, 5e-15);
    assert!(
        result
            .stage
            .final_inventory
            .iter()
            .all(|(_, amount)| amount.value() >= 0.0)
    );
}

#[test]
fn zero_reagent_or_rates_need_no_numerical_steps_but_still_check_horizon() {
    for (k, acid) in [(0.005, 0.0), (0.0, 20.0)] {
        let network = network(vec![copper(k), dolomite(k)]);
        let result = run(&network, acid, 10.0);
        assert_eq!(result.stage.final_inventory, inventory(&network, acid));
        assert_eq!(result.accepted_steps, 0);
        assert_eq!(result.rejected_steps, 0);
        assert!(
            result
                .stage
                .reactions
                .iter()
                .all(|r| r.extent.value() == 0.0)
        );
        assert!(matches!(
            network.advance(
                &inventory(&network, acid),
                Volume::new(1.0),
                Time::new(101.0),
                Temperature::new(298.15),
                &conditions(&network),
                budget()
            ),
            Err(BatchNetworkError::Source { .. })
        ));
    }
}

#[test]
fn work_and_source_domain_refusals_are_atomic_and_use_total_elapsed_time() {
    let network = network(vec![copper(0.005), dolomite(0.005)]);
    let initial = inventory(&network, 40.0);
    let saved = initial.clone();
    let mut short_budget = budget();
    short_budget.maximum_attempts = 1;
    short_budget.maximum_step = Time::new(0.05);
    assert!(matches!(
        network.advance(
            &initial,
            Volume::new(1.0),
            Time::new(10.0),
            Temperature::new(298.15),
            &conditions(&network),
            short_budget
        ),
        Err(BatchNetworkError::BudgetExhausted {
            attempts: 1,
            accepted_steps: 1,
            ..
        })
    ));
    assert!(matches!(
        network.advance(
            &initial,
            Volume::new(1.0),
            Time::new(101.0),
            Temperature::new(298.15),
            &conditions(&network),
            budget()
        ),
        Err(BatchNetworkError::Source { .. })
    ));
    assert_eq!(initial, saved);

    let base = copper(0.005);
    let mut narrow = card(0.005);
    narrow.validity = narrow.validity.with_quantity(
        REACTANT_A_CONCENTRATION_AXIS,
        QuantitySpec::dimensional(AmountConcentration::DIMS),
        9.5,
        10.0,
    );
    let restricted = BilinearBatchReaction::from_card(
        base.reaction().clone(),
        base.reactants()[0].clone(),
        base.reactants()[1].clone(),
        narrow,
    )
    .unwrap();
    let network =
        super::BilinearBatchNetwork::new(vec![restricted], SpeciesId::new("H+(aq)").unwrap())
            .unwrap();
    let initial = inventory(&network, 20.0);
    let saved = initial.clone();
    let mut limited = budget();
    limited.maximum_attempts = 24;
    assert!(matches!(
        network.advance(
            &initial,
            Volume::new(1.0),
            Time::new(10.0),
            Temperature::new(298.15),
            &conditions(&network),
            limited
        ),
        Err(BatchNetworkError::BudgetExhausted { .. }) | Err(BatchNetworkError::Source { .. })
    ));
    assert_eq!(initial, saved);
}

#[test]
fn missing_conditions_owned_axis_spoofing_and_unsupported_topologies_refuse() {
    let network = network(vec![copper(0.005), dolomite(0.005)]);
    let initial = inventory(&network, 40.0);
    assert!(matches!(
        network.advance(
            &initial,
            Volume::new(1.0),
            Time::new(1.0),
            Temperature::new(298.15),
            &conditions(&network)[..1],
            budget()
        ),
        Err(BatchNetworkError::Invalid(_))
    ));
    let mut queries = conditions(&network);
    queries[0].1 = QueryPoint::new().with(ELAPSED_TIME_AXIS, 0.0).unwrap();
    assert!(matches!(
        network.advance(
            &initial,
            Volume::new(1.0),
            Time::new(1.0),
            Temperature::new(298.15),
            &queries,
            budget()
        ),
        Err(BatchNetworkError::Invalid(_))
    ));
    assert!(
        BilinearBatchNetwork::new(
            vec![copper(0.005), dolomite(0.005)],
            SpeciesId::new("CuO(s)").unwrap()
        )
        .is_err()
    );
    assert!(
        BilinearBatchNetwork::new(
            vec![copper(0.005), copper(0.006)],
            SpeciesId::new("H+(aq)").unwrap()
        )
        .is_err()
    );
    let mut invalid = budget();
    invalid.maximum_attempts = 0;
    assert!(matches!(
        network.advance(
            &initial,
            Volume::new(1.0),
            Time::new(1.0),
            Temperature::new(298.15),
            &conditions(&network),
            invalid
        ),
        Err(BatchNetworkError::Invalid(_))
    ));
}

#[test]
fn positive_subnormal_rates_do_not_disappear_in_the_rk_weighted_mean() {
    let smallest = f64::from_bits(1);
    let duration = 1e308;
    let base = copper(smallest);
    let mut source = card(smallest);
    source.validity = source.validity.with_quantity(
        ELAPSED_TIME_AXIS,
        QuantitySpec::dimensional(Time::DIMS),
        0.0,
        duration,
    );
    let model = BilinearBatchReaction::from_card(
        base.reaction().clone(),
        base.reactants()[0].clone(),
        base.reactants()[1].clone(),
        source,
    )
    .unwrap();
    let network = network(vec![model]);
    let initial = network
        .reaction_network()
        .species()
        .iter()
        .map(|id| {
            let value = if matches!(id.as_str(), "CuO(s)" | "H+(aq)") {
                1.0
            } else {
                0.0
            };
            (id.clone(), Amount::new(value))
        })
        .collect::<Vec<_>>();
    let controls = BatchIntegrationBudget {
        maximum_step: Time::new(duration),
        absolute_extent_error: Amount::new(1e-14),
        relative_extent_error: Dimensionless::new(0.0),
        maximum_attempts: 1,
    };
    let result = network
        .advance(
            &initial,
            Volume::new(1.0),
            Time::new(duration),
            Temperature::new(298.15),
            &conditions(&network),
            controls,
        )
        .unwrap();
    // All four derivatives round to the least positive binary64 number. Its
    // weighted mean is still that number, and its finite-time integral is
    // large enough to change the actual mineral and proton inventories.
    close(
        result.stage.reactions[0].extent.value(),
        smallest * duration,
        1e-30,
    );
    assert!(result.stage.reactions[0].extent.value() > 0.0);
    assert!(amount(&result, "CuO(s)") < 1.0);
    assert!(amount(&result, "H+(aq)") < 1.0);
}

#[test]
fn nonzero_error_difference_cannot_underflow_to_a_zero_error_claim() {
    assert!(matches!(
        local_error_estimate(f64::from_bits(1), 0.0),
        Err(BatchNetworkError::Arithmetic(
            "positive local-error estimate range loss"
        ))
    ));
}
