//! Continuous-flow CIL limits, retained carbon and open-stream balances.

use super::*;
use crate::cyanide_leach::{
    cil::{CarbonRole, IonPairTransfer, tests as cil_fixture},
    flow::tests as flow_fixture,
    tests as chemistry_fixture,
};

fn active_tank(rate: f64, carbon_rate: f64) -> CarbonInLeachCstr {
    let charge = chemistry_fixture::charge(1.0);
    let tank = CyanidationCstr::new(
        flow_fixture::model(rate, false),
        charge,
        QueryPoint::new(),
        chemistry_fixture::options(),
    )
    .unwrap();
    CarbonInLeachCstr::from_tank(
        tank,
        IonPairTransfer::sodium().unwrap(),
        vec![
            cil_fixture::population(
                "product-carbon",
                CarbonRole::ProductCarbon,
                1.0,
                0.0,
                0.6,
                0.1,
                carbon_rate,
                10.0,
            ),
            cil_fixture::population(
                "preg-robber",
                CarbonRole::PregRobbing,
                1.0,
                0.0,
                0.3,
                0.2,
                0.25 * carbon_rate,
                10.0,
            ),
        ],
    )
    .unwrap()
}

fn balanced(answer: &FlowCilAdvance) {
    for balance in [
        &answer.interval_conservation,
        &answer.cumulative_conservation,
    ] {
        flow_fixture::close(balance.gold.value(), 0.0, 4e-10);
        flow_fixture::close(balance.counterion.value(), 0.0, 4e-10);
        flow_fixture::close(balance.dissolved_oxygen.value(), 0.0, 4e-10);
        flow_fixture::close(balance.maximum_site_residual.value(), 0.0, 4e-10);
        let ligand = answer
            .state
            .solution
            .competition
            .as_ref()
            .map_or(0.0, |c| c.ligand_balance_residual.value());
        flow_fixture::close(balance.cyanide.value(), ligand, 4e-9);
    }
}

#[test]
fn zero_carbon_rate_matches_original_continuous_leach() {
    let charge = chemistry_fixture::charge(1.0);
    let feed = flow_fixture::feed(&charge);
    let base = CyanidationCstr::new(
        flow_fixture::model(4.0, false),
        charge,
        QueryPoint::new(),
        chemistry_fixture::options(),
    )
    .unwrap();
    let mut expected = base.clone();
    let expected = expected
        .advance(
            Time::new(0.4),
            VolumetricFlowRate::new(0.2),
            &feed,
            flow_fixture::numerics(),
            || false,
        )
        .unwrap();
    let mut actual = CarbonInLeachCstr::from_tank(
        base,
        IonPairTransfer::sodium().unwrap(),
        vec![cil_fixture::population(
            "inert-carbon",
            CarbonRole::ProductCarbon,
            1.0,
            0.0,
            0.6,
            0.1,
            0.0,
            10.0,
        )],
    )
    .unwrap();
    let answer = actual
        .advance(
            Time::new(0.4),
            VolumetricFlowRate::new(0.2),
            &feed,
            flow_fixture::numerics(),
            || false,
        )
        .unwrap();
    for (actual, expected) in [
        (
            answer.state.solution.inventory.accessible_gold.value(),
            expected.state.inventory.accessible_gold.value(),
        ),
        (
            answer.state.solution.inventory.cyanide_ligand.value(),
            expected.state.inventory.cyanide_ligand.value(),
        ),
        (
            answer.state.solution.inventory.dissolved_oxygen.value(),
            expected.state.inventory.dissolved_oxygen.value(),
        ),
        (
            answer.state.solution.inventory.gold_complex.value(),
            expected.state.inventory.gold_complex.value(),
        ),
        (
            answer.interval_gold_dissolved.value(),
            expected.interval_gold_dissolved.value(),
        ),
        (
            answer.outlet.gold_complex.value(),
            expected.outlet.gold_complex.value(),
        ),
    ] {
        flow_fixture::close(actual, expected, 4e-8);
    }
    assert_eq!(answer.interval_product_capture.value(), 0.0);
    assert_eq!(answer.interval_preg_robbing_capture.value(), 0.0);
    balanced(&answer);
}

#[test]
fn flowing_solution_and_retained_carbon_capture_gold_simultaneously() {
    let charge = chemistry_fixture::charge(1.0);
    let feed = flow_fixture::feed(&charge);
    let mut tank = active_tank(4.0, 1.4);
    let answer = tank
        .advance(
            Time::new(0.6),
            VolumetricFlowRate::new(0.2),
            &feed,
            flow_fixture::numerics(),
            || false,
        )
        .unwrap();
    assert!(answer.interval_gold_dissolved.value() > 0.0);
    assert!(answer.interval_product_capture.value() > 0.0);
    assert!(answer.interval_preg_robbing_capture.value() > 0.0);
    assert!(
        answer.interval_product_capture.value()
            > answer.interval_preg_robbing_capture.value()
    );
    assert!(answer.outlet.gold_complex.value() > 0.0);
    assert!(answer.state.solution.inventory.gold_complex.value() > 0.0);
    assert_eq!(answer.state.populations.len(), 2);
    assert_eq!(tank.models().len(), 2);
    balanced(&answer);

    let original_charge = tank.tank().initial_charge().clone();
    let continued = tank.advance(
        Time::new(0.1), VolumetricFlowRate::new(0.1), &feed,
        flow_fixture::numerics(), || false,
    ).unwrap();
    balanced(&continued);
    assert_eq!(tank.tank().initial_charge(), &original_charge);
    flow_fixture::close(continued.previous_elapsed.value(), 0.6, 1e-15);
    flow_fixture::close(continued.cumulative_inlet.accessible_gold.value(), 0.13, 2e-15);
    flow_fixture::close(
        continued.cumulative_outlet.gold_complex.value(),
        answer.cumulative_outlet.gold_complex.value() + continued.outlet.gold_complex.value(),
        2e-15,
    );
    let stored: f64 = continued.state.populations.iter().map(|p| p.adsorbed.value()).sum();
    let mobile = &continued.state.solution.inventory;
    let outgoing = &continued.cumulative_outlet;
    let incoming = &continued.cumulative_inlet;
    flow_fixture::close(
        mobile.accessible_gold.value() + mobile.gold_complex.value() + stored
            + outgoing.accessible_gold.value() + outgoing.gold_complex.value(),
        original_charge.accessible_gold.value() + original_charge.gold_complex.value()
            + incoming.accessible_gold.value() + incoming.gold_complex.value(),
        4e-10,
    );
    flow_fixture::close(
        mobile.cyanide_ligand.value() + 2.0 * (mobile.gold_complex.value() + stored)
            + outgoing.cyanide_ligand.value() + 2.0 * outgoing.gold_complex.value(),
        original_charge.uncomplexed_cyanide.value() + 2.0 * original_charge.gold_complex.value()
            + incoming.cyanide_ligand.value() + 2.0 * incoming.gold_complex.value(),
        4e-10,
    );
    let sodium = |inventory: &CyanidationInventory| inventory.spectators.iter()
        .find(|ion| ion.species.as_str() == "Na[+](aq)").unwrap().amount.value();
    let initial_sodium = original_charge.spectators.iter()
        .find(|ion| ion.species.as_str() == "Na[+](aq)").unwrap().amount.value();
    flow_fixture::close(
        sodium(mobile) + stored + sodium(outgoing),
        initial_sodium + sodium(incoming),
        4e-10,
    );
}

#[test]
fn stopped_liquid_flow_reduces_to_retained_carbon_batch_physics() {
    let charge = chemistry_fixture::charge(1.0);
    let feed = flow_fixture::feed(&charge);
    let mut tank = active_tank(4.0, 1.4);
    let answer = tank
        .advance(
            Time::new(0.25),
            VolumetricFlowRate::new(0.0),
            &feed,
            flow_fixture::numerics(),
            || false,
        )
        .unwrap();
    assert!(answer.interval_gold_dissolved.value() > 0.0);
    assert!(answer.interval_product_capture.value() > 0.0);
    assert_eq!(answer.inlet.accessible_gold.value(), 0.0);
    assert_eq!(answer.outlet.accessible_gold.value(), 0.0);
    assert_eq!(answer.outlet.gold_complex.value(), 0.0);
    balanced(&answer);
}

#[test]
fn counterion_feed_guard_is_atomic() {
    let charge = chemistry_fixture::charge(1.0);
    let mut feed = flow_fixture::feed(&charge);
    let mut tank = active_tank(4.0, 1.4);
    let before = format!("{:?}", tank.state());
    feed.spectators
        .iter_mut()
        .find(|ion| ion.species.as_str() == "Na[+](aq)")
        .unwrap()
        .concentration = AmountConcentration::new(0.1);
    assert!(matches!(
        tank.advance(
            Time::new(0.2),
            VolumetricFlowRate::new(0.2),
            &feed,
            flow_fixture::numerics(),
            || false,
        ),
        Err(FlowCilError::Invalid(
            "counterion feed must cover all incoming accessible and complexed Au"
        ))
    ));
    assert_eq!(format!("{:?}", tank.state()), before);
}

#[test]
fn zero_carbon_rate_still_admits_original_sources_and_preserves_typed_refusal() {
    use crate::aqueous_adsorption::AqueousAdsorptionError;
    let charge = chemistry_fixture::charge(1.0);
    let feed = flow_fixture::feed(&charge);
    let base = CyanidationCstr::new(
        flow_fixture::model(4.0, false), charge, QueryPoint::new(),
        chemistry_fixture::options(),
    ).unwrap();
    let mut tank = CarbonInLeachCstr::from_tank(
        base,
        IonPairTransfer::sodium().unwrap(),
        vec![cil_fixture::population(
            "inert-source-qualified-carbon", CarbonRole::ProductCarbon,
            1.0, 0.0, 0.6, 0.1, 0.0, 1e-8,
        )],
    ).unwrap();
    let before = format!("{tank:?}");
    assert!(matches!(
        tank.advance(
            Time::new(0.1), VolumetricFlowRate::new(0.2), &feed,
            flow_fixture::numerics(), || false,
        ),
        Err(FlowCilError::Adsorption(AdsorptionFlowError::Source(
            AqueousAdsorptionError::OutsideValidity { .. }
        )))
    ));
    assert_eq!(format!("{tank:?}"), before);
}
