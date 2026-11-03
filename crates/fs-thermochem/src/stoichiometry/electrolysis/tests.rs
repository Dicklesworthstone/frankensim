use super::*;
use crate::stoichiometry::{FARADAY_C_PER_MOL, minerals::MineralReaction};

fn id(name: &str) -> SpeciesId {
    SpeciesId::new(name).unwrap()
}
fn model() -> ElectrolysisReaction {
    ElectrolysisReaction::new(
        MineralReaction::CopperElectrowinning.build().unwrap(),
        id("Cu(s)"),
        MolarMass::new(0.063546),
        2,
    )
    .unwrap()
}
fn inventory(copper: f64, water: f64) -> Vec<(SpeciesId, Amount)> {
    [
        ("Cu[2+](aq)", copper),
        ("H2O(l)", water),
        ("Cu(s)", 0.0),
        ("O2(g)", 0.0),
        ("H+(aq)", 0.0),
    ]
    .into_iter()
    .map(|(name, value)| (id(name), Amount::new(value)))
    .collect()
}
fn value(result: &ElectrolysisAdvance, species: &str) -> f64 {
    result
        .modeled_final_inventory
        .iter()
        .find(|(name, _)| name.as_str() == species)
        .unwrap()
        .1
        .value()
}
fn close(actual: f64, expected: f64) {
    assert!((actual - expected).abs() <= expected.abs().max(1.0) * 2e-14);
}

#[test]
fn copper_net_equation_preserves_stocks_atoms_charge_and_electron_basis() {
    let model = model();
    let initial = inventory(3.0, 5.0);
    let preserved = initial.clone();
    let result = model
        .advance_charge(
            &initial,
            ElectricCharge::new(2.0 * FARADAY_C_PER_MOL),
            Dimensionless::new(1.0),
        )
        .unwrap();
    close(result.extent.value(), 0.5);
    close(result.product_yield.amount().value(), 1.0);
    close(result.product_yield.mass().value(), 0.063546);
    close(value(&result, "Cu[2+](aq)"), 2.0);
    close(value(&result, "H2O(l)"), 4.0);
    close(value(&result, "Cu(s)"), 1.0);
    close(value(&result, "O2(g)"), 0.5);
    close(value(&result, "H+(aq)"), 2.0);
    close(value(&result, "Cu[2+](aq)") + value(&result, "Cu(s)"), 3.0);
    close(
        2.0 * value(&result, "H2O(l)") + value(&result, "H+(aq)"),
        10.0,
    );
    close(
        value(&result, "H2O(l)") + 2.0 * value(&result, "O2(g)"),
        5.0,
    );
    close(
        2.0 * value(&result, "Cu[2+](aq)") + value(&result, "H+(aq)"),
        6.0,
    );
    assert_eq!(model.electrons_per_product(), 2);
    assert_eq!(model.electrons_per_equation(), 4);
    assert_eq!(result.model.reaction().id(), model.reaction().id());
    assert_eq!(result.model.reaction().terms(), model.reaction().terms());
    assert_eq!(result.model.product(), &id("Cu(s)"));
    assert_eq!(result.model.molar_mass(), MolarMass::new(0.063546));
    assert_eq!(result.unallocated_charge.value(), 0.0);
    assert_eq!(initial, preserved);
    let mut reversed = initial.clone();
    reversed.reverse();
    assert_eq!(
        model
            .advance_charge(&reversed, result.passed_charge, Dimensionless::new(1.0))
            .unwrap()
            .modeled_final_inventory,
        result.modeled_final_inventory,
    );
}

#[test]
fn partial_efficiency_reports_unallocated_charge_without_inventing_side_products() {
    let result = model()
        .advance_charge(
            &inventory(3.0, 5.0),
            ElectricCharge::new(4.0 * FARADAY_C_PER_MOL),
            Dimensionless::new(0.75),
        )
        .unwrap();
    close(result.extent.value(), 0.75);
    close(result.product_yield.amount().value(), 1.5);
    close(result.productive_charge.value(), 3.0 * FARADAY_C_PER_MOL);
    close(result.unallocated_charge.value(), FARADAY_C_PER_MOL);
    close(
        result.productive_charge.value() + result.unallocated_charge.value(),
        result.passed_charge.value(),
    );
    // Only the productive pathway's water, acid and oxygen are accounted for.
    close(value(&result, "H2O(l)"), 3.5);
    close(value(&result, "H+(aq)"), 3.0);
    close(value(&result, "O2(g)"), 0.75);
    assert_eq!(result.modeled_final_inventory.len(), 5);
}

#[test]
fn constant_current_known_yield_and_zero_drive_keep_explicit_inventory_semantics() {
    let model = model();
    let initial = inventory(100.0, 1000.0);
    let result = model
        .advance_current(
            &initial,
            Current::new(1000.0),
            Time::new(3600.0),
            Dimensionless::new(1.0),
        )
        .unwrap();
    close(result.passed_charge.value(), 3_600_000.0);
    close(result.product_yield.mass().value(), 1.185494183238305);
    let direct = model
        .advance_charge(
            &initial,
            ElectricCharge::new(3_600_000.0),
            Dimensionless::new(1.0),
        )
        .unwrap();
    assert_eq!(
        result.modeled_final_inventory,
        direct.modeled_final_inventory
    );
    for (charge, efficiency) in [(0.0, 1.0), (100.0, 0.0)] {
        let result = model
            .advance_charge(
                &initial,
                ElectricCharge::new(charge),
                Dimensionless::new(efficiency),
            )
            .unwrap();
        assert_eq!(result.extent.value(), 0.0);
        assert_eq!(result.product_yield.mass().value(), 0.0);
        let mut ordered = initial.clone();
        ordered.sort_by(|a, b| a.0.cmp(&b.0));
        assert_eq!(result.modeled_final_inventory, ordered);
        assert_eq!(result.unallocated_charge.value(), charge);
    }
}

#[test]
fn charge_cannot_overdraw_copper_or_water_and_failures_are_atomic() {
    let model = model();
    let initial = inventory(2.0, 2.0);
    let preserved = initial.clone();
    let charge = 4.0 * FARADAY_C_PER_MOL;
    let result = model
        .advance_charge(
            &initial,
            ElectricCharge::new(charge),
            Dimensionless::new(1.0),
        )
        .unwrap();
    assert_eq!(value(&result, "Cu[2+](aq)"), 0.0);
    assert_eq!(value(&result, "H2O(l)"), 0.0);
    for (stock, charge) in [
        (initial.clone(), charge.next_up()),
        (inventory(10.0, 0.5), charge),
    ] {
        assert!(matches!(
            model.advance_charge(&stock, ElectricCharge::new(charge), Dimensionless::new(1.0)),
            Err(ReactionBalanceError::ExtentExceedsInventory { .. })
        ));
    }
    assert_eq!(initial, preserved);
    let mut duplicate = initial.clone();
    duplicate[1] = duplicate[0].clone();
    assert!(
        model
            .advance_charge(
                &duplicate,
                ElectricCharge::new(0.0),
                Dimensionless::new(1.0)
            )
            .is_err()
    );
    assert!(
        model
            .advance_charge(
                &initial[..4],
                ElectricCharge::new(0.0),
                Dimensionless::new(1.0)
            )
            .is_err()
    );
}

#[test]
fn invalid_associations_drives_and_unrepresentable_products_refuse() {
    let reaction = MineralReaction::CopperElectrowinning.build().unwrap();
    for product in ["Cu[2+](aq)", "unknown"] {
        assert!(
            ElectrolysisReaction::new(reaction.clone(), id(product), MolarMass::new(0.063546), 2)
                .is_err()
        );
    }
    assert!(
        ElectrolysisReaction::new(reaction.clone(), id("Cu(s)"), MolarMass::new(0.063546), 0)
            .is_err()
    );
    for mass in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            ElectrolysisReaction::new(reaction.clone(), id("Cu(s)"), MolarMass::new(mass), 2)
                .is_err()
        );
    }
    let model = model();
    let initial = inventory(100.0, 1000.0);
    for bad in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            model
                .advance_charge(&initial, ElectricCharge::new(bad), Dimensionless::new(1.0))
                .is_err()
        );
        assert!(
            model
                .advance_charge(&initial, ElectricCharge::new(1.0), Dimensionless::new(bad))
                .is_err()
        );
        assert!(
            model
                .advance_current(
                    &initial,
                    Current::new(bad),
                    Time::new(1.0),
                    Dimensionless::new(1.0)
                )
                .is_err()
        );
        assert!(
            model
                .advance_current(
                    &initial,
                    Current::new(1.0),
                    Time::new(bad),
                    Dimensionless::new(1.0)
                )
                .is_err()
        );
    }
    assert!(
        model
            .advance_charge(&initial, ElectricCharge::new(1.0), Dimensionless::new(1.01))
            .is_err()
    );
    for (current, time) in [(f64::MAX, 2.0), (f64::from_bits(1), f64::from_bits(1))] {
        assert!(
            model
                .advance_current(
                    &initial,
                    Current::new(current),
                    Time::new(time),
                    Dimensionless::new(1.0)
                )
                .is_err()
        );
    }
    assert!(
        model
            .advance_charge(
                &initial,
                ElectricCharge::new(f64::from_bits(1)),
                Dimensionless::new(1.0)
            )
            .is_err()
    );
}
