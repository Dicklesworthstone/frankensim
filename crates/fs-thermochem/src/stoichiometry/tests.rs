//! G0 exact conservation and physical amount identities; no kinetic/field claim.
use super::minerals::MineralReaction;
use super::*;
use crate::ElementId;

fn id(value: &str) -> SpeciesId {
    SpeciesId::new(value).unwrap()
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-15 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}
fn amount(reaction: &ConservedReaction, extent: f64, species: &str) -> f64 {
    reaction
        .amounts(Amount::new(extent))
        .unwrap()
        .into_iter()
        .find(|a| a.term.species == id(species))
        .unwrap()
        .amount
        .value()
}
fn water(coefficients: [i128; 3], charges: [i128; 3]) -> Result<ConservedReaction> {
    let species = vec![id("H2"), id("O2"), id("H2O")];
    let a = ElementalMatrix::new(
        vec![ElementId::new("H").unwrap(), ElementId::new("O").unwrap()],
        species.clone(),
        vec![vec![2, 0, 2], vec![0, 2, 1]],
    )?;
    let reaction = ReactionId::new("water").unwrap();
    let n = StoichiometricMatrix::new(
        species.clone(),
        vec![reaction.clone()],
        coefficients.into_iter().map(|n| vec![n]).collect(),
    )?;
    let z = ChargeVector::new(species, charges.to_vec())?;
    ConservedReaction::new(&a, &n, &z, reaction)
}

#[test]
fn all_reference_reactions_conserve_and_rebuild_identically() {
    use MineralReaction::*;
    for recipe in [
        LimeHydration,
        CalciteAcidNeutralization,
        DolomiteAcidNeutralization,
        PyriteOxidationFerrous,
        PyriteOxidationFerricHydroxide,
    ] {
        let first = recipe.build().unwrap();
        let second = recipe.build().unwrap();
        assert_eq!(
            first.conservation().identity(),
            second.conservation().identity()
        );
        assert!(
            first
                .terms()
                .windows(2)
                .all(|w| w[0].species < w[1].species)
        );
        assert!(
            first
                .amounts(Amount::new(0.0))
                .unwrap()
                .iter()
                .all(|a| a.amount.value() == 0.0)
        );
    }
}

#[test]
fn pyrite_endpoint_changes_oxygen_demand_and_acid_generation() {
    let ferrous = MineralReaction::PyriteOxidationFerrous.build().unwrap();
    let ferric = MineralReaction::PyriteOxidationFerricHydroxide
        .build()
        .unwrap();
    for (reaction, extent, oxygen, acid, water) in [
        (&ferrous, 0.5, 3.5, 2.0, 1.0),
        (&ferric, 0.25, 3.75, 4.0, 3.5),
    ] {
        close(amount(reaction, extent, "FeS2(s)"), 1.0);
        close(amount(reaction, extent, "O2(g)"), oxygen);
        close(amount(reaction, extent, "H+(aq)"), acid);
        close(amount(reaction, extent, "H2O(l)"), water);
        close(amount(reaction, extent, "SO4[2-](aq)"), 2.0);
    }
    // An explicit caller mass basis of 0.120 kg/mol: 0.240 kg reacted is 2 mol FeS2.
    close(
        ferric
            .extent_from_reacted_mass(&id("FeS2(s)"), Mass::new(0.240), MolarMass::new(0.120))
            .unwrap()
            .value(),
        0.5,
    );
}

#[test]
fn carbonate_neutralization_and_lime_amounts() {
    let calcite = MineralReaction::CalciteAcidNeutralization.build().unwrap();
    let dolomite = MineralReaction::DolomiteAcidNeutralization.build().unwrap();
    close(amount(&calcite, 1.0, "H+(aq)"), 2.0);
    close(amount(&calcite, 1.0, "CO2(aq)"), 1.0);
    close(amount(&dolomite, 1.0, "H+(aq)"), 4.0);
    close(amount(&dolomite, 1.0, "CO2(aq)"), 2.0);
    close(amount(&dolomite, 1.0, "Mg[2+](aq)"), 1.0);
    let lime = MineralReaction::LimeHydration.build().unwrap();
    close(amount(&lime, 3.0, "CaO(s)"), 3.0);
    close(amount(&lime, 3.0, "H2O(l)"), 3.0);
    close(amount(&lime, 3.0, "Ca(OH)2(s)"), 3.0);
}

#[test]
fn independent_element_charge_and_degenerate_equation_refusals() {
    assert!(matches!(
        water([-1, -1, 2], [0, 0, 0]),
        Err(ReactionBalanceError::Chemistry(
            ChemistryError::ElementImbalance { .. }
        ))
    ));
    assert!(matches!(
        water([-2, -1, 2], [0, 1, 0]),
        Err(ReactionBalanceError::Chemistry(
            ChemistryError::ChargeImbalance { .. }
        ))
    ));
    assert!(water([0, 0, 0], [0, 0, 0]).is_err());
    let oversized = (1_i128 << 53) + 1;
    assert!(water([-2 * oversized, -oversized, 2 * oversized], [0, 0, 0]).is_err());
    assert!(water([-2, -1, 2], [0, 0, 0]).is_ok());
}

#[test]
fn limiting_inventory_uses_each_reactant_once_and_refuses_missing_values() {
    let reaction = water([-2, -1, 2], [0, 0, 0]).unwrap();
    let inventory = |h, o| vec![(id("H2"), Amount::new(h)), (id("O2"), Amount::new(o))];
    close(
        reaction
            .limiting_extent(&inventory(6.0, 2.0))
            .unwrap()
            .value(),
        2.0,
    );
    close(
        reaction
            .limiting_extent(&inventory(2.0, 6.0))
            .unwrap()
            .value(),
        1.0,
    );
    close(
        reaction
            .limiting_extent(&inventory(0.0, 2.0))
            .unwrap()
            .value(),
        0.0,
    );
    assert!(reaction.limiting_extent(&[]).is_err());
    assert!(
        reaction
            .limiting_extent(&[(id("H2"), Amount::new(1.0))])
            .is_err()
    );
    assert!(
        reaction
            .limiting_extent(&[(id("H2"), Amount::new(1.0)), (id("H2"), Amount::new(1.0))])
            .is_err()
    );
    assert!(
        reaction
            .limiting_extent(&[(id("H2"), Amount::new(1.0)), (id("H2O"), Amount::new(1.0))])
            .is_err()
    );
    for bad in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(reaction.limiting_extent(&inventory(bad, 2.0)).is_err());
        assert!(reaction.amounts(Amount::new(bad)).is_err());
        assert!(
            reaction
                .extent_from_reacted_mass(&id("H2"), Mass::new(bad), MolarMass::new(0.002))
                .is_err()
        );
    }
    assert!(
        reaction
            .extent_from_reacted_mass(&id("H2O"), Mass::new(1.0), MolarMass::new(0.018))
            .is_err()
    );
}

#[test]
fn overflow_and_positive_underflow_refuse_instead_of_zero_amounts() {
    let reaction = water([-2, -1, 2], [0, 0, 0]).unwrap();
    assert!(reaction.amounts(Amount::new(f64::MAX)).is_err());
    assert!(
        reaction
            .extent_from_reacted_mass(
                &id("H2"),
                Mass::new(f64::MAX),
                MolarMass::new(f64::MIN_POSITIVE)
            )
            .is_err()
    );
    assert!(
        reaction
            .limiting_extent(&[
                (id("H2"), Amount::new(f64::from_bits(1))),
                (id("O2"), Amount::new(1.0))
            ])
            .is_err()
    );
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            reaction
                .extent_from_reacted_mass(&id("H2"), Mass::new(1.0), MolarMass::new(bad))
                .is_err()
        );
    }
}

#[test]
fn faraday_copper_yield_uses_electron_count_and_explicit_efficiency() {
    let molar_mass = MolarMass::new(0.063546);
    let full = faraday_yield(
        ElectricCharge::new(2.0 * FARADAY_C_PER_MOL),
        molar_mass,
        2,
        1.0,
    )
    .unwrap();
    close(full.amount().value(), 1.0);
    close(full.mass().value(), 0.063546);
    let partial = faraday_yield(
        ElectricCharge::new(2.0 * FARADAY_C_PER_MOL),
        molar_mass,
        2,
        0.85,
    )
    .unwrap();
    close(partial.mass().value(), 0.063546 * 0.85);
    close(
        faraday_yield(ElectricCharge::new(3_600_000.0), molar_mass, 2, 1.0)
            .unwrap()
            .mass()
            .value(),
        1.185494183238305,
    );
    assert_eq!(
        faraday_yield(ElectricCharge::new(1.0), molar_mass, 2, 0.0)
            .unwrap()
            .mass()
            .value(),
        0.0
    );
    assert_eq!(
        faraday_yield(ElectricCharge::new(0.0), molar_mass, 2, 1.0)
            .unwrap()
            .amount()
            .value(),
        0.0
    );
}

#[test]
fn faraday_invalid_charge_efficiency_mass_and_range_refuse() {
    for charge in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            faraday_yield(
                ElectricCharge::new(charge),
                MolarMass::new(0.063546),
                2,
                1.0
            )
            .is_err()
        );
    }
    for efficiency in [-0.1, 1.01, f64::NAN, f64::INFINITY] {
        assert!(
            faraday_yield(
                ElectricCharge::new(1.0),
                MolarMass::new(0.063546),
                2,
                efficiency
            )
            .is_err()
        );
    }
    for mass in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(faraday_yield(ElectricCharge::new(1.0), MolarMass::new(mass), 2, 1.0).is_err());
    }
    assert!(faraday_yield(ElectricCharge::new(1.0), MolarMass::new(0.063546), 0, 1.0).is_err());
    assert!(
        faraday_yield(
            ElectricCharge::new(f64::MAX),
            MolarMass::new(f64::MAX),
            2,
            1.0
        )
        .is_err()
    );
    assert!(
        faraday_yield(
            ElectricCharge::new(f64::from_bits(1)),
            MolarMass::new(0.063546),
            2,
            1.0
        )
        .is_err()
    );
}

#[test]
fn enthalpy_extent_sign_zero_and_range_are_preserved() {
    let heat = ReactionEnthalpy {
        molar: MolarEnergyQuantityV1::new(-65_000.0),
        terms: vec![],
    };
    assert_eq!(
        heat.at_extent(Amount::new(2.0)).unwrap().value(),
        -130_000.0
    );
    assert_eq!(heat.at_extent(Amount::new(0.0)).unwrap().value(), 0.0);
    for bad in [-1.0, f64::NAN, f64::INFINITY, f64::MAX] {
        assert!(heat.at_extent(Amount::new(bad)).is_err());
    }
}
