use super::*;
use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{LawId, LawParameter, Provenance};
use fs_thermochem::stoichiometry::minerals::MineralReaction;

fn id(s: &str) -> SpeciesId {
    SpeciesId::new(s).unwrap()
}
fn card(k: f64) -> ConstitutiveModelCard {
    ConstitutiveModelCard {
        law: LawId(LOCAL_BILINEAR_LAW.into()),
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
            "local-bilinear-cstr-fixture",
            b"authored rates; no mineral calibration",
        )],
        provenance: Provenance {
            source: "synthetic local CSTR kinetics".into(),
            license: "project fixture".into(),
            artifact: None,
        },
    }
}
fn copper(k: f64) -> LocalBilinearReaction {
    LocalBilinearReaction::from_card(
        MineralReaction::CupricOxideAcidLeach.build().unwrap(),
        id("CuO(s)"),
        id("H+(aq)"),
        card(k),
    )
    .unwrap()
}
fn calcite(k: f64) -> LocalBilinearReaction {
    LocalBilinearReaction::from_card(
        MineralReaction::CalciteBicarbonateNeutralization
            .build()
            .unwrap(),
        id("CaCO3(s)"),
        id("H+(aq)"),
        card(k),
    )
    .unwrap()
}
fn tank(models: Vec<LocalBilinearReaction>) -> ReactiveCstr {
    ReactiveCstr::new(Volume::new(1.0), models, id("H+(aq)")).unwrap()
}
fn inventory(tank: &ReactiveCstr, values: &[(&str, f64)]) -> Vec<(SpeciesId, Amount)> {
    tank.reaction_network()
        .species()
        .iter()
        .map(|s| {
            let n = values
                .iter()
                .find(|(name, _)| *name == s.as_str())
                .map_or(0.0, |(_, n)| *n);
            (s.clone(), Amount::new(n))
        })
        .collect()
}
fn feed(tank: &ReactiveCstr, values: &[(&str, f64)]) -> Vec<(SpeciesId, AmountConcentration)> {
    inventory(tank, values)
        .into_iter()
        .map(|(s, n)| (s, AmountConcentration::new(n.value())))
        .collect()
}
fn conditions(tank: &ReactiveCstr) -> Vec<(ReactionId, QueryPoint)> {
    tank.models()
        .iter()
        .map(|m| (m.reaction().id().clone(), QueryPoint::new()))
        .collect()
}
fn budget() -> FlowIntegrationBudget {
    FlowIntegrationBudget {
        maximum_step: Time::new(1.0),
        absolute_amount_error: Amount::new(1e-8),
        relative_amount_error: Dimensionless::new(1e-7),
        maximum_attempts: 8192,
    }
}
fn run(
    tank: &ReactiveCstr,
    initial: &[(SpeciesId, Amount)],
    inlet: &[(SpeciesId, AmountConcentration)],
    q: f64,
    t: f64,
    numerics: FlowIntegrationBudget,
) -> ReactionFlowAdvance {
    tank.advance(
        initial,
        inlet,
        VolumetricFlowRate::new(q),
        Time::new(t),
        Temperature::new(298.15),
        &conditions(tank),
        numerics,
    )
    .unwrap()
}
fn amount(values: &[(SpeciesId, Amount)], species: &str) -> f64 {
    values
        .iter()
        .find(|(s, _)| s.as_str() == species)
        .unwrap()
        .1
        .value()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a} != {b} (tolerance {tolerance})"
    );
}
fn balances(result: &ReactionFlowAdvance, tolerance: f64) {
    for (species, residual) in &result.balance_residuals {
        assert!(
            residual.value().abs() <= tolerance,
            "{} residual {}",
            species.as_str(),
            residual.value()
        );
    }
    for (i, (s, n)) in result.final_inventory.iter().enumerate() {
        assert!(n.value() >= 0.0);
        assert_eq!(s, &result.inlet_amounts[i].0);
        close(
            result.initial_inventory[i].1.value()
                + result.inlet_amounts[i].1.value()
                + result.reactions.formed[i].1.value(),
            n.value() + result.outlet_amounts[i].1.value() + result.reactions.consumed[i].1.value(),
            tolerance,
        );
    }
}

#[test]
fn stopped_flow_matches_independent_bilinear_closed_form() {
    let model = copper(0.005);
    let tank = tank(vec![model.clone()]);
    let initial = inventory(&tank, &[("CuO(s)", 10.0), ("H+(aq)", 20.0)]);
    let result = run(&tank, &initial, &feed(&tank, &[]), 0.0, 10.0, budget());
    // H=2 CuO throughout; dn_CuO/dt=-.01 n_CuO², n(10)=10/(1+.01*10*10)=5.
    close(amount(&result.final_inventory, "CuO(s)"), 5.0, 2e-6);
    close(amount(&result.final_inventory, "H+(aq)"), 10.0, 4e-6);
    close(amount(&result.final_inventory, "Cu[2+](aq)"), 5.0, 2e-6);
    close(result.reactions.reactions[0].extent.value(), 5.0, 2e-6);
    assert!(result.outlet_amounts.iter().all(|(_, n)| n.value() == 0.0));
    assert_eq!(result.model().models()[0].source(), model.source());
    assert_eq!(
        result.model().models()[0].reaction().terms(),
        model.reaction().terms()
    );
    assert!(result.accepted_steps > 1);
    assert!(result.estimated_extent_errors[0].1.value() <= 1.01e-6);
    assert!(amount(&result.estimated_inventory_errors, "CuO(s)") <= 1.01e-6);
    balances(&result, 1e-12);
}

#[test]
fn zero_rate_matches_exact_nonreactive_flow_and_independent_outlet_integral() {
    let tank = tank(vec![copper(0.0)]);
    let initial = inventory(
        &tank,
        &[
            ("CuO(s)", 10.0),
            ("H+(aq)", 4.0),
            ("Cu[2+](aq)", 2.0),
            ("H2O(l)", 100.0),
        ],
    );
    let inlet = feed(
        &tank,
        &[
            ("CuO(s)", 3.0),
            ("H+(aq)", 8.0),
            ("Cu[2+](aq)", 5.0),
            ("H2O(l)", 1000.0),
        ],
    );
    let result = run(&tank, &initial, &inlet, 0.2, 3.0, budget());
    let survival = (-0.6_f64).exp();
    for (i, (_, n)) in result.final_inventory.iter().enumerate() {
        let initial = initial[i].1.value();
        let target = inlet[i].1.value();
        close(n.value(), target + (initial - target) * survival, 2e-13);
        close(
            result.outlet_amounts[i].1.value(),
            initial * (1.0 - survival) + target * (0.6 - 1.0 + survival),
            4e-13,
        );
    }
    assert_eq!(result.accepted_steps, 1);
    assert_eq!(result.reactions.reactions[0].extent.value(), 0.0);
    assert!(
        result
            .estimated_outlet_errors
            .iter()
            .all(|(_, n)| n.value() == 0.0)
    );
    balances(&result, 3e-13);
}

#[test]
fn competing_paths_share_finite_acid_with_different_stoichiometric_demands() {
    let tank = tank(vec![copper(0.01), calcite(0.002)]);
    let initial = inventory(&tank, &[("CuO(s)", 10.0), ("CaCO3(s)", 10.0)]);
    let stopped = run(&tank, &initial, &feed(&tank, &[]), 0.1, 2.0, budget());
    assert!(
        stopped
            .reactions
            .reactions
            .iter()
            .all(|r| r.extent.value() == 0.0)
    );
    let result = run(
        &tank,
        &stopped.final_inventory,
        &feed(&tank, &[("H+(aq)", 8.0)]),
        0.1,
        10.0,
        budget(),
    );
    assert!(
        result
            .reactions
            .reactions
            .iter()
            .all(|r| r.extent.value() > 0.0)
    );
    let copper_extent = result
        .reactions
        .reactions
        .iter()
        .find(|r| {
            r.amounts
                .iter()
                .any(|a| a.term.species.as_str() == "CuO(s)")
        })
        .unwrap()
        .extent
        .value();
    let calcite_extent = result
        .reactions
        .reactions
        .iter()
        .find(|r| {
            r.amounts
                .iter()
                .any(|a| a.term.species.as_str() == "CaCO3(s)")
        })
        .unwrap()
        .extent
        .value();
    let acid_used = 2.0 * copper_extent + calcite_extent;
    assert!(acid_used > 0.0 && acid_used < 8.0);
    close(
        acid_used
            + amount(&result.final_inventory, "H+(aq)")
            + amount(&result.outlet_amounts, "H+(aq)"),
        8.0,
        3e-12,
    );
    assert!(amount(&result.outlet_amounts, "Cu[2+](aq)") > 0.0);
    assert!(amount(&result.final_inventory, "Cu[2+](aq)") > 0.0);
    // Stopping the pump does not stop chemistry in the retained wet inventory.
    let hold = run(
        &tank,
        &result.final_inventory,
        &feed(&tank, &[]),
        0.0,
        3.0,
        budget(),
    );
    assert!(amount(&hold.final_inventory, "H+(aq)") < amount(&result.final_inventory, "H+(aq)"));
    assert!(hold.outlet_amounts.iter().all(|(_, n)| n.value() == 0.0));
    balances(&result, 5e-12);
    balances(&hold, 5e-12);
    // Independently check integrated element and charge totals across the tank.
    for component in 0..6 {
        let counts = |name: &str| -> [f64; 6] {
            match name {
                "CuO(s)" => [1., 0., 0., 1., 0., 0.],
                "Cu[2+](aq)" => [1., 0., 0., 0., 0., 2.],
                "CaCO3(s)" => [0., 1., 1., 3., 0., 0.],
                "Ca[2+](aq)" => [0., 1., 0., 0., 0., 2.],
                "HCO3[-](aq)" => [0., 0., 1., 3., 1., -1.],
                "H+(aq)" => [0., 0., 0., 0., 1., 1.],
                "H2O(l)" => [0., 0., 0., 1., 2., 0.],
                _ => panic!("unexpected species {name}"),
            }
        };
        let total = |values: &[(SpeciesId, Amount)]| {
            values
                .iter()
                .map(|(s, n)| counts(s.as_str())[component] * n.value())
                .sum::<f64>()
        };
        close(
            total(&result.initial_inventory) + total(&result.inlet_amounts),
            total(&result.final_inventory) + total(&result.outlet_amounts),
            2e-11,
        );
    }
}

#[test]
fn stationary_balance_and_time_partition_refine_at_second_order() {
    let tank = tank(vec![copper(0.005)]);
    // .2(10-x)=.01x², H=2x, Cu2+=water=10-x at the exact unsplit steady state.
    let x = (-0.2_f64 + (0.12_f64).sqrt()) / 0.02;
    let initial = inventory(
        &tank,
        &[
            ("CuO(s)", x),
            ("H+(aq)", 2.0 * x),
            ("Cu[2+](aq)", 10.0 - x),
            ("H2O(l)", 10.0 - x),
        ],
    );
    let inlet = feed(&tank, &[("CuO(s)", 10.0), ("H+(aq)", 20.0)]);
    let mut coarse_budget = budget();
    coarse_budget.absolute_amount_error = Amount::new(1.0);
    coarse_budget.relative_amount_error = Dimensionless::new(0.0);
    coarse_budget.maximum_step = Time::new(0.4);
    let coarse = run(&tank, &initial, &inlet, 0.2, 8.0, coarse_budget);
    let mut fine_budget = coarse_budget;
    fine_budget.maximum_step = Time::new(0.2);
    let fine = run(&tank, &initial, &inlet, 0.2, 8.0, fine_budget);
    let coarse_error = (amount(&coarse.final_inventory, "CuO(s)") - x).abs();
    let fine_error = (amount(&fine.final_inventory, "CuO(s)") - x).abs();
    assert!(
        coarse_error > 3.5 * fine_error && coarse_error < 4.5 * fine_error,
        "{coarse_error} / {fine_error}"
    );
    assert!(fine_error < 1e-3);
    // Two separately committed physical intervals converge to the same continuation.
    let tight = run(&tank, &initial, &inlet, 0.2, 8.0, budget());
    let first = run(&tank, &initial, &inlet, 0.2, 4.0, budget());
    let second = run(&tank, &first.final_inventory, &inlet, 0.2, 4.0, budget());
    for (i, (_, n)) in tight.final_inventory.iter().enumerate() {
        close(n.value(), second.final_inventory[i].1.value(), 5e-6);
        close(
            tight.outlet_amounts[i].1.value(),
            first.outlet_amounts[i].1.value() + second.outlet_amounts[i].1.value(),
            5e-6,
        );
    }
    close(amount(&tight.final_inventory, "CuO(s)"), x, 5e-6);
    balances(&tight, 2e-11);
}

#[test]
fn model_inventory_feed_and_condition_permutations_preserve_results() {
    let a = tank(vec![copper(0.01), calcite(0.002)]);
    let b = tank(vec![calcite(0.002), copper(0.01)]);
    let initial = inventory(&a, &[("CuO(s)", 10.0), ("CaCO3(s)", 10.0), ("H+(aq)", 3.0)]);
    let inlet = feed(&a, &[("H+(aq)", 8.0)]);
    let first = run(&a, &initial, &inlet, 0.1, 2.0, budget());
    let mut reordered_initial = initial.clone();
    reordered_initial.reverse();
    let mut reordered_feed = inlet.clone();
    reordered_feed.reverse();
    let mut reordered_conditions = conditions(&b);
    reordered_conditions.reverse();
    let second = b
        .advance(
            &reordered_initial,
            &reordered_feed,
            VolumetricFlowRate::new(0.1),
            Time::new(2.0),
            Temperature::new(298.15),
            &reordered_conditions,
            budget(),
        )
        .unwrap();
    assert_eq!(first.final_inventory, second.final_inventory);
    assert_eq!(first.outlet_amounts, second.outlet_amounts);
    assert_eq!(
        first.estimated_extent_errors,
        second.estimated_extent_errors
    );
    assert_eq!(first.accepted_steps, second.accepted_steps);
}

#[test]
fn local_law_refuses_batch_age_residence_dimensions_and_incomplete_axes() {
    let build = |source| {
        LocalBilinearReaction::from_card(
            MineralReaction::CupricOxideAcidLeach.build().unwrap(),
            id("CuO(s)"),
            id("H+(aq)"),
            source,
        )
    };
    let mut batch = card(0.01);
    batch.law = LawId("bilinear-batch-reaction".into());
    assert!(build(batch).is_err());
    for axis in ["elapsed-time", "mean-residence-time"] {
        let mut source = card(0.01);
        source.validity =
            source
                .validity
                .with_quantity(axis, QuantitySpec::dimensional(Time::DIMS), 0.0, 10.0);
        assert!(build(source).is_err());
    }
    let mut wrong = card(0.01);
    wrong.parameters.get_mut(RATE_PARAMETER).unwrap().dims = Time::DIMS;
    assert!(build(wrong).is_err());
    let tank = tank(vec![copper(0.01)]);
    let initial = inventory(&tank, &[("CuO(s)", 10.0), ("H+(aq)", 3.0)]);
    let mut inlet = feed(&tank, &[]);
    inlet.pop();
    assert!(
        tank.advance(
            &initial,
            &inlet,
            VolumetricFlowRate::new(0.1),
            Time::new(2.0),
            Temperature::new(298.15),
            &conditions(&tank),
            budget()
        )
        .is_err()
    );
}

#[test]
fn exhausted_work_and_late_source_departure_return_no_partial_state() {
    let mut source = card(0.01);
    source.validity = source.validity.with_quantity(
        REACTANT_B_CONCENTRATION_AXIS,
        QuantitySpec::dimensional(AmountConcentration::DIMS),
        0.0,
        1.0,
    );
    let limited = LocalBilinearReaction::from_card(
        MineralReaction::CupricOxideAcidLeach.build().unwrap(),
        id("CuO(s)"),
        id("H+(aq)"),
        source,
    )
    .unwrap();
    let tank = tank(vec![limited]);
    let initial = inventory(&tank, &[("CuO(s)", 10.0)]);
    let saved = initial.clone();
    let mut short = budget();
    short.maximum_attempts = 32;
    assert!(
        tank.advance(
            &initial,
            &feed(&tank, &[("H+(aq)", 50.0)]),
            VolumetricFlowRate::new(0.1),
            Time::new(10.0),
            Temperature::new(298.15),
            &conditions(&tank),
            short
        )
        .is_err()
    );
    assert_eq!(initial, saved);
    let mut tiny = budget();
    tiny.maximum_attempts = 1;
    assert!(matches!(
        tank.advance(
            &initial,
            &feed(&tank, &[]),
            VolumetricFlowRate::new(0.1),
            Time::new(10.0),
            Temperature::new(298.15),
            &conditions(&tank),
            tiny
        ),
        Err(ReactionFlowError::BudgetExhausted { .. })
    ));
    assert_eq!(initial, saved);
}
