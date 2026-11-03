use super::*;
use crate::stoichiometry::minerals::MineralReaction;
use crate::{ChargeVector, ElementId, ElementalMatrix, StoichiometricMatrix};

fn species(name: &str) -> SpeciesId {
    SpeciesId::new(name).unwrap()
}
fn reaction(name: &str) -> ReactionId {
    ReactionId::new(name).unwrap()
}
fn initial(network: &ReactionNetwork, values: &[(&str, f64)]) -> Vec<(SpeciesId, Amount)> {
    network
        .species()
        .iter()
        .map(|id| {
            (
                id.clone(),
                Amount::new(
                    values
                        .iter()
                        .find(|(name, _)| *name == id.as_str())
                        .map_or(0.0, |(_, amount)| *amount),
                ),
            )
        })
        .collect()
}
fn extents(network: &ReactionNetwork, values: &[(&str, f64)]) -> Vec<(ReactionId, Amount)> {
    network
        .reactions()
        .iter()
        .map(|r| {
            (
                r.id().clone(),
                Amount::new(
                    values
                        .iter()
                        .find(|(name, _)| *name == r.id().as_str())
                        .map_or(0.0, |(_, amount)| *amount),
                ),
            )
        })
        .collect()
}
fn value(values: &[(SpeciesId, Amount)], name: &str) -> f64 {
    values
        .iter()
        .find(|(id, _)| id.as_str() == name)
        .unwrap()
        .1
        .value()
}
fn acid_network() -> ReactionNetwork {
    ReactionNetwork::new(vec![
        MineralReaction::MalachiteAcidLeach.build().unwrap(),
        MineralReaction::CalciteAcidNeutralization.build().unwrap(),
    ])
    .unwrap()
}

// Independent synthetic isomers: conservation stays admitted by fs-qty.
fn isomer(
    id: &str,
    from: &str,
    to: &str,
    atoms: i128,
    charge: i128,
    zero_axis: bool,
) -> ConservedReaction {
    let ids = vec![species(from), species(to)];
    let mut elements = vec![ElementId::new("H").unwrap()];
    let mut counts = vec![vec![atoms, atoms]];
    if zero_axis {
        elements.push(ElementId::new("O").unwrap());
        counts.push(vec![0, 0]);
    }
    let elemental = ElementalMatrix::new(elements, ids.clone(), counts).unwrap();
    let stoichiometric =
        StoichiometricMatrix::new(ids.clone(), vec![reaction(id)], vec![vec![-1], vec![1]])
            .unwrap();
    let charges = ChargeVector::new(ids, vec![charge, charge]).unwrap();
    ConservedReaction::new(&elemental, &stoichiometric, &charges, reaction(id)).unwrap()
}

#[test]
fn shared_acid_competes_and_an_admitted_stage_conserves_atoms_and_charge() {
    let network = acid_network();
    let inventory = initial(
        &network,
        &[("CaCO3(s)", 1.0), ("Cu2CO3(OH)2(s)", 1.0), ("H+(aq)", 4.0)],
    );
    let untouched = inventory.clone();
    let excessive = extents(
        &network,
        &[("calcite-acid-co2", 1.0), ("malachite-acid-leach", 1.0)],
    );
    assert!(matches!(network.advance_inventory(&inventory, &excessive),
        Err(ReactionNetworkError::InsufficientInventory { species: id, available, required })
        if id == species("H+(aq)") && available.value() == 4.0 && required.value() == 6.0));
    let stage = network
        .advance_inventory(
            &inventory,
            &extents(
                &network,
                &[("calcite-acid-co2", 1.0), ("malachite-acid-leach", 0.5)],
            ),
        )
        .unwrap();
    for (name, amount) in [
        ("CaCO3(s)", 0.0),
        ("Cu2CO3(OH)2(s)", 0.5),
        ("H+(aq)", 0.0),
        ("Ca[2+](aq)", 1.0),
        ("Cu[2+](aq)", 1.0),
        ("CO2(aq)", 1.5),
        ("H2O(l)", 2.5),
    ] {
        assert_eq!(value(&stage.final_inventory, name), amount);
    }
    assert_eq!(value(&stage.consumed, "H+(aq)"), 4.0);
    assert_eq!(value(&stage.formed, "H2O(l)"), 2.5);
    assert_eq!(stage.reactions[0].reaction, reaction("calcite-acid-co2"));
    assert_eq!(stage.reactions[1].extent.value(), 0.5);
    // Independent elemental/charge accounts for these seven physical species.
    let totals = |pool: &[(SpeciesId, Amount)]| {
        let ca = value(pool, "CaCO3(s)");
        let cu = value(pool, "Cu2CO3(OH)2(s)");
        let h = value(pool, "H+(aq)");
        let ca_ion = value(pool, "Ca[2+](aq)");
        let cu_ion = value(pool, "Cu[2+](aq)");
        let co2 = value(pool, "CO2(aq)");
        let water = value(pool, "H2O(l)");
        [
            ca + ca_ion,
            2.0 * cu + cu_ion,
            ca + cu + co2,
            3.0 * ca + 5.0 * cu + 2.0 * co2 + water,
            2.0 * cu + h + 2.0 * water,
            h + 2.0 * ca_ion + 2.0 * cu_ion,
        ]
    };
    assert_eq!(totals(&inventory), totals(&stage.final_inventory));
    assert_eq!(inventory, untouched);
}

#[test]
fn pyrite_acid_can_fund_only_a_later_explicit_neutralization_stage() {
    let network = ReactionNetwork::new(vec![
        MineralReaction::PyriteOxidationFerrous.build().unwrap(),
        MineralReaction::CalciteAcidNeutralization.build().unwrap(),
    ])
    .unwrap();
    let inventory = initial(
        &network,
        &[
            ("FeS2(s)", 2.0),
            ("O2(g)", 7.0),
            ("H2O(l)", 2.0),
            ("CaCO3(s)", 2.0),
        ],
    );
    let together = extents(
        &network,
        &[("pyrite-oxidation-ferrous", 1.0), ("calcite-acid-co2", 2.0)],
    );
    assert!(matches!(network.advance_inventory(&inventory, &together),
        Err(ReactionNetworkError::InsufficientInventory { species: id, .. }) if id == species("H+(aq)")));
    let oxidation = network
        .advance_inventory(
            &inventory,
            &extents(&network, &[("pyrite-oxidation-ferrous", 1.0)]),
        )
        .unwrap();
    assert_eq!(value(&oxidation.final_inventory, "H+(aq)"), 4.0);
    assert_eq!(value(&oxidation.final_inventory, "O2(g)"), 0.0);
    assert_eq!(value(&oxidation.final_inventory, "H2O(l)"), 0.0);
    let neutralization = network
        .advance_inventory(
            &oxidation.final_inventory,
            &extents(&network, &[("calcite-acid-co2", 2.0)]),
        )
        .unwrap();
    assert_eq!(value(&neutralization.final_inventory, "H+(aq)"), 0.0);
    assert_eq!(value(&neutralization.final_inventory, "H2O(l)"), 2.0);
    assert_eq!(value(&neutralization.final_inventory, "Ca[2+](aq)"), 2.0);
    assert_eq!(value(&neutralization.final_inventory, "SO4[2-](aq)"), 4.0);
}

#[test]
fn species_definitions_match_nonzero_atoms_and_charge_not_incidental_zero_axes() {
    let first = isomer("first", "A", "B", 1, 0, false);
    for (atoms, charge) in [(2, 0), (1, 1)] {
        assert!(
            matches!(ReactionNetwork::new(vec![first.clone(), isomer("second", "B", "C", atoms, charge, false)]),
            Err(ReactionNetworkError::ConflictingSpecies { species: id }) if id == species("B"))
        );
    }
    let network =
        ReactionNetwork::new(vec![isomer("second", "B", "C", 1, 0, true), first.clone()]).unwrap();
    assert_eq!(
        network.species(),
        &[species("A"), species("B"), species("C")]
    );
    assert_eq!(network.reactions()[0].id(), &reaction("first"));
    assert!(ReactionNetwork::new(vec![first.clone(), first]).is_err());
    assert!(ReactionNetwork::new(Vec::new()).is_err());
}

#[test]
fn explicit_turnover_with_initial_stock_is_valid_but_zero_stock_cannot_fund_a_cycle() {
    let network = ReactionNetwork::new(vec![
        isomer("forward", "A", "B", 1, 0, false),
        isomer("reverse", "B", "A", 1, 0, false),
    ])
    .unwrap();
    let both = extents(&network, &[("forward", 1.0), ("reverse", 1.0)]);
    assert!(matches!(
        network.advance_inventory(&initial(&network, &[]), &both),
        Err(ReactionNetworkError::InsufficientInventory { .. })
    ));
    let inventory = initial(&network, &[("A", 1.0), ("B", 1.0)]);
    let stage = network.advance_inventory(&inventory, &both).unwrap();
    assert_eq!(stage.final_inventory, inventory);
    assert_eq!(value(&stage.consumed, "A"), 1.0);
    assert_eq!(value(&stage.formed, "A"), 1.0);
}

#[test]
fn axes_are_complete_unique_and_finite_and_permutations_preserve_results() {
    let network = acid_network();
    let inventory = initial(
        &network,
        &[("CaCO3(s)", 1.0), ("Cu2CO3(OH)2(s)", 1.0), ("H+(aq)", 6.0)],
    );
    let requested = extents(
        &network,
        &[("calcite-acid-co2", 0.5), ("malachite-acid-leach", 0.25)],
    );
    let expected = network.advance_inventory(&inventory, &requested).unwrap();
    let (mut reversed_inventory, mut reversed_requested) = (inventory.clone(), requested.clone());
    reversed_inventory.reverse();
    reversed_requested.reverse();
    let reordered = network
        .advance_inventory(&reversed_inventory, &reversed_requested)
        .unwrap();
    assert_eq!(reordered.final_inventory, expected.final_inventory);
    assert_eq!(reordered.consumed, expected.consumed);
    assert_eq!(reordered.formed, expected.formed);
    assert_eq!(
        network
            .advance_inventory(&inventory, &extents(&network, &[]))
            .unwrap()
            .final_inventory,
        inventory
    );
    assert!(
        network
            .advance_inventory(&inventory[..inventory.len() - 1], &requested)
            .is_err()
    );
    assert!(
        network
            .advance_inventory(&inventory, &requested[..1])
            .is_err()
    );
    for duplicate in [true, false] {
        let mut invalid = inventory.clone();
        invalid[0].0 = if duplicate {
            invalid[1].0.clone()
        } else {
            species("unknown")
        };
        assert!(network.advance_inventory(&invalid, &requested).is_err());
        let mut invalid = requested.clone();
        invalid[0].0 = if duplicate {
            invalid[1].0.clone()
        } else {
            reaction("unknown")
        };
        assert!(network.advance_inventory(&inventory, &invalid).is_err());
    }
    for bad in [-1.0, f64::NAN, f64::INFINITY] {
        let mut invalid = inventory.clone();
        invalid[0].1 = Amount::new(bad);
        assert!(network.advance_inventory(&invalid, &requested).is_err());
        let mut invalid = requested.clone();
        invalid[0].1 = Amount::new(bad);
        assert!(network.advance_inventory(&inventory, &invalid).is_err());
    }
}

#[test]
fn rounded_product_and_shared_demand_cannot_hide_stock_overdraw() {
    let network = ReactionNetwork::new(vec![
        MineralReaction::PyriteOxidationFerricHydroxide
            .build()
            .unwrap(),
    ])
    .unwrap();
    let inventory = initial(
        &network,
        &[("FeS2(s)", 1.0), ("O2(g)", 0.1), ("H2O(l)", 1.0)],
    );
    let nearest = 0.1_f64 / 15.0;
    assert_eq!(nearest * 15.0, 0.1);
    assert!(nearest.mul_add(15.0, -0.1) > 0.0);
    let id = network.reactions()[0].id().clone();
    assert!(
        matches!(network.advance_inventory(&inventory, &[(id.clone(), Amount::new(nearest))]),
        Err(ReactionNetworkError::InsufficientInventory { species: id, .. }) if id == species("O2(g)"))
    );
    let feasible = nearest.next_down();
    let stage = network
        .advance_inventory(&inventory, &[(id, Amount::new(feasible))])
        .unwrap();
    assert_eq!(
        value(&stage.final_inventory, "O2(g)"),
        (-15.0_f64).mul_add(feasible, 0.1)
    );

    let shared = ReactionNetwork::new(vec![
        isomer("one", "A", "B", 1, 0, false),
        isomer("two", "A", "C", 1, 0, false),
    ])
    .unwrap();
    let inventory = initial(&shared, &[("A", 1.0)]);
    let requested = extents(&shared, &[("one", 1.0), ("two", f64::from_bits(1))]);
    assert!(matches!(shared.advance_inventory(&inventory, &requested),
        Err(ReactionNetworkError::InsufficientInventory { available, required, .. })
        if available.value() == 1.0 && required.value() == 1.0));
}

#[test]
fn nonzero_resolution_loss_and_gross_overflow_refuse_but_minimal_amounts_survive() {
    let network = ReactionNetwork::new(vec![isomer("one", "A", "B", 1, 0, false)]).unwrap();
    assert!(matches!(
        network.advance_inventory(
            &initial(&network, &[("A", 1.0), ("B", 1e300)]),
            &extents(&network, &[("one", 1.0)])
        ),
        Err(ReactionNetworkError::Arithmetic(_))
    ));
    let tiny = f64::from_bits(1);
    let stage = network
        .advance_inventory(
            &initial(&network, &[("A", tiny)]),
            &extents(&network, &[("one", tiny)]),
        )
        .unwrap();
    assert_eq!(value(&stage.final_inventory, "A"), 0.0);
    assert_eq!(value(&stage.final_inventory, "B"), tiny);
    let network = ReactionNetwork::new(vec![
        isomer("one", "A", "C", 1, 0, false),
        isomer("two", "B", "C", 1, 0, false),
    ])
    .unwrap();
    assert!(matches!(
        network.advance_inventory(
            &initial(&network, &[("A", 1e308), ("B", 1e308)]),
            &extents(&network, &[("one", 1e308), ("two", 1e308)])
        ),
        Err(ReactionNetworkError::Arithmetic(_))
    ));
}
