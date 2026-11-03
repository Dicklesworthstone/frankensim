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
        CalciteCarbonicAcidNeutralization,
        DolomiteAcidNeutralization,
        CalciteBicarbonateNeutralization,
        DolomiteBicarbonateNeutralization,
        HydratedLimeNeutralization,
        MalachiteAcidLeach,
        MalachiteCarbonicAcidLeach,
        GoldCyanidation,
        GoldCyanidationDissolvedOxygen,
        HydrogenCyanideDissociation,
        HydrogenSulfateDissociation,
        WaterDissociation,
        PyriteOxidationFerrous,
        FerrousOxidationFerric,
        FerricPyriteOxidation,
        PyriteOxidationFerricHydroxide,
        MagnesiteCalcination,
        MagnetiteOxidationToHematite,
        CuprousOxideOxidation,
        CupricOxideCarbonMonoxideReduction,
        CuprousOxideCarbonMonoxideReduction,
        PyriteRoastingToHematite,
        SulfurDioxideOxidation,
        SulfurTrioxideHydration,
        PyriteToHematiteAndSulfuricAcid,
        CopperElectrowinning,
        FerricChalcopyriteLeach,
        FerricChalcociteLeach,
        FerricCovelliteLeach,
        FerricBorniteLeach,
        SphaleriteRoastingToZincite,
        ZinciteAcidLeach,
        ZincElectrowinning,
        FerricSphaleriteLeach,
        CopperCementationWithZinc,
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
fn zinc_roast_leach_electrowinning_closes_intermediates_and_respects_stage_stocks() {
    let roast = MineralReaction::SphaleriteRoastingToZincite
        .build()
        .unwrap();
    let leach = MineralReaction::ZinciteAcidLeach.build().unwrap();
    let electrowinning = MineralReaction::ZincElectrowinning.build().unwrap();
    // Two mol sphalerite through the three separately declared stages yields
    // 2 ZnS + 2 O2 -> 2 Zn + 2 SO2. Acid/water cancel only across stages.
    let mut net = std::collections::BTreeMap::new();
    for (reaction, extent) in [(&roast, 1), (&leach, 2), (&electrowinning, 1)] {
        for term in reaction.terms() {
            *net.entry(term.species.as_str()).or_insert(0_i128) += extent * term.coefficient;
        }
    }
    net.retain(|_, coefficient| *coefficient != 0);
    assert_eq!(
        net,
        std::collections::BTreeMap::from([
            ("ZnS(s)", -2),
            ("O2(g)", -2),
            ("Zn(s)", 2),
            ("SO2(g)", 2),
        ])
    );
    close(amount(&roast, 0.5, "O2(g)"), 1.5);
    close(amount(&leach, 1.0, "H+(aq)"), 2.0);
    close(amount(&electrowinning, 0.5, "H+(aq)"), 2.0);
    let initial = [
        (id("ZnO(s)"), Amount::new(3.0)),
        (id("H+(aq)"), Amount::new(4.0)),
        (id("Zn[2+](aq)"), Amount::new(0.0)),
        (id("H2O(l)"), Amount::new(0.0)),
    ];
    let available = leach.limiting_extent(&initial[..2]).unwrap();
    close(available.value(), 2.0);
    let final_inventory = leach.advance_inventory(&initial, available).unwrap();
    for ((species, _), expected) in initial.iter().zip([1.0, 0.0, 2.0, 2.0]) {
        close(
            final_inventory
                .iter()
                .find(|(id, _)| id == species)
                .unwrap()
                .1
                .value(),
            expected,
        );
    }
    assert!(
        leach
            .advance_inventory(&initial, Amount::new(2.01))
            .is_err()
    );
    let ferric = MineralReaction::FerricSphaleriteLeach.build().unwrap();
    close(amount(&ferric, 1.0, "Fe[3+](aq)"), 2.0);
    close(amount(&ferric, 1.0, "Zn[2+](aq)"), 1.0);
    close(amount(&ferric, 1.0, "Fe[2+](aq)"), 2.0);
    close(amount(&ferric, 1.0, "S(s)"), 1.0);
    let cementation = MineralReaction::CopperCementationWithZinc.build().unwrap();
    for species in ["Cu[2+](aq)", "Zn(s)", "Cu(s)", "Zn[2+](aq)"] {
        close(amount(&cementation, 1.5, species), 1.5);
    }
}

#[test]
fn ferric_copper_sulfides_keep_mineral_iron_and_oxidant_stocks_separate() {
    use MineralReaction::*;
    for (recipe, mineral, ferric, copper, ferrous, sulfur) in [
        (FerricChalcopyriteLeach, "CuFeS2(s)", 4.0, 1.0, 5.0, 2.0),
        (FerricChalcociteLeach, "Cu2S(s)", 4.0, 2.0, 4.0, 1.0),
        (FerricCovelliteLeach, "CuS(s)", 2.0, 1.0, 2.0, 1.0),
        (FerricBorniteLeach, "Cu5FeS4(s)", 12.0, 5.0, 13.0, 4.0),
    ] {
        let reaction = recipe.build().unwrap();
        close(amount(&reaction, 1.0, mineral), 1.0);
        close(amount(&reaction, 1.0, "Fe[3+](aq)"), ferric);
        close(amount(&reaction, 1.0, "Cu[2+](aq)"), copper);
        close(amount(&reaction, 1.0, "Fe[2+](aq)"), ferrous);
        close(amount(&reaction, 1.0, "S(s)"), sulfur);
        // Enough mineral for three extents but ferric for only two. No
        // regenerated ferric is borrowed from a hypothetical later stage.
        let initial = [
            (id(mineral), Amount::new(3.0)),
            (id("Fe[3+](aq)"), Amount::new(2.0 * ferric)),
            (id("Cu[2+](aq)"), Amount::new(0.0)),
            (id("Fe[2+](aq)"), Amount::new(0.0)),
            (id("S(s)"), Amount::new(0.0)),
        ];
        let extent = reaction.limiting_extent(&initial[..2]).unwrap();
        close(extent.value(), 2.0);
        let advanced = reaction.advance_inventory(&initial, extent).unwrap();
        for ((species, _), expected) in
            initial
                .iter()
                .zip([1.0, 0.0, 2.0 * copper, 2.0 * ferrous, 2.0 * sulfur])
        {
            close(
                advanced
                    .iter()
                    .find(|(s, _)| s == species)
                    .unwrap()
                    .1
                    .value(),
                expected,
            );
        }
        assert!(
            reaction
                .advance_inventory(&initial, Amount::new(2.01))
                .is_err()
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
fn ferric_pyrite_cycle_cancels_intermediate_iron_and_regeneration_acid() {
    let attack = MineralReaction::FerricPyriteOxidation.build().unwrap();
    let regeneration = MineralReaction::FerrousOxidationFerric.build().unwrap();
    let net = MineralReaction::PyriteOxidationFerrous.build().unwrap();
    // Two ferric attacks + seven integer regeneration equations equal the
    // direct two-pyrite endpoint. Do not count 16 H+ per pyrite as net acid.
    let mut combined = std::collections::BTreeMap::new();
    for (reaction, multiple) in [(&attack, 2), (&regeneration, 7), (&net, -1)] {
        for term in reaction.terms() {
            *combined.entry(term.species.clone()).or_insert(0_i128) += multiple * term.coefficient;
        }
    }
    assert!(combined.values().all(|coefficient| *coefficient == 0));
    close(amount(&attack, 1.0, "Fe[3+](aq)"), 14.0);
    close(amount(&attack, 1.0, "Fe[2+](aq)"), 15.0);
    close(amount(&attack, 1.0, "H+(aq)"), 16.0);
    close(amount(&regeneration, 3.5, "H+(aq)"), 14.0);
    close(amount(&regeneration, 3.5, "O2(g)"), 3.5);
}

#[test]
fn ferric_pyrite_leach_stops_at_available_oxidant() {
    let reaction = MineralReaction::FerricPyriteOxidation.build().unwrap();
    let initial = [
        (id("FeS2(s)"), Amount::new(3.0)),
        (id("Fe[3+](aq)"), Amount::new(28.0)),
        (id("H2O(l)"), Amount::new(20.0)),
        (id("Fe[2+](aq)"), Amount::new(0.0)),
        (id("SO4[2-](aq)"), Amount::new(0.0)),
        (id("H+(aq)"), Amount::new(0.0)),
    ];
    let extent = reaction.limiting_extent(&initial[..3]).unwrap();
    close(extent.value(), 2.0);
    let final_inventory = reaction.advance_inventory(&initial, extent).unwrap();
    for ((species, _), expected) in initial.iter().zip([1.0, 0.0, 4.0, 30.0, 4.0, 32.0]) {
        close(
            final_inventory
                .iter()
                .find(|(id, _)| id == species)
                .unwrap()
                .1
                .value(),
            expected,
        );
    }
    assert!(
        reaction
            .advance_inventory(&initial, Amount::new(2.01))
            .is_err()
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
fn carbonic_acid_endpoints_retain_distinct_species_and_water_balance() {
    let calcite = MineralReaction::CalciteCarbonicAcidNeutralization
        .build()
        .unwrap();
    let malachite = MineralReaction::MalachiteCarbonicAcidLeach.build().unwrap();
    close(amount(&calcite, 1.0, "H+(aq)"), 2.0);
    close(amount(&calcite, 1.0, "H2CO3(aq)"), 1.0);
    assert!(
        !calcite
            .terms()
            .iter()
            .any(|term| term.species == id("H2O(l)"))
    );
    close(amount(&malachite, 1.0, "H+(aq)"), 4.0);
    close(amount(&malachite, 1.0, "Cu[2+](aq)"), 2.0);
    close(amount(&malachite, 1.0, "H2CO3(aq)"), 1.0);
    close(amount(&malachite, 1.0, "H2O(l)"), 2.0);
    for reaction in [calcite, malachite] {
        assert!(
            !reaction
                .terms()
                .iter()
                .any(|term| term.species == id("CO2(aq)"))
        );
    }
}

#[test]
fn calcination_and_oxide_oxidation_keep_gas_endpoints_and_written_extent() {
    let calcination = MineralReaction::MagnesiteCalcination.build().unwrap();
    close(amount(&calcination, 2.0, "MgCO3(s)"), 2.0);
    close(amount(&calcination, 2.0, "MgO(s)"), 2.0);
    close(amount(&calcination, 2.0, "CO2(g)"), 2.0);
    let oxidation = MineralReaction::MagnetiteOxidationToHematite
        .build()
        .unwrap();
    close(amount(&oxidation, 0.25, "Fe3O4(s)"), 1.0);
    close(amount(&oxidation, 0.25, "O2(g)"), 0.25);
    close(amount(&oxidation, 0.25, "Fe2O3(s)"), 1.5);
    // Oxygen limits this authored inventory: 10 mol magnetite needs 2.5 mol O2.
    let extent = oxidation
        .limiting_extent(&[
            (id("Fe3O4(s)"), Amount::new(10.0)),
            (id("O2(g)"), Amount::new(2.0)),
        ])
        .unwrap();
    close(extent.value(), 2.0);
    close(amount(&oxidation, extent.value(), "Fe3O4(s)"), 8.0);
}

#[test]
fn roasting_and_acid_endpoints_preserve_the_sulfur_and_oxygen_basis() {
    let roasting = MineralReaction::PyriteRoastingToHematite.build().unwrap();
    // One mole pyrite is a quarter mole of the written four-pyrite equation.
    close(amount(&roasting, 0.25, "FeS2(s)"), 1.0);
    close(amount(&roasting, 0.25, "O2(g)"), 2.75);
    close(amount(&roasting, 0.25, "Fe2O3(s)"), 0.5);
    close(amount(&roasting, 0.25, "SO2(g)"), 2.0);
    assert!(
        roasting
            .terms()
            .iter()
            .all(|term| !term.species.as_str().contains("aq"))
    );
    let available = roasting
        .limiting_extent(&[
            (id("FeS2(s)"), Amount::new(4.0)),
            (id("O2(g)"), Amount::new(5.5)),
        ])
        .unwrap();
    close(available.value(), 0.5);
    close(amount(&roasting, available.value(), "FeS2(s)"), 2.0);

    let oxidation = MineralReaction::SulfurDioxideOxidation.build().unwrap();
    close(amount(&oxidation, 1.0, "SO2(g)"), 2.0);
    close(amount(&oxidation, 1.0, "O2(g)"), 1.0);
    close(amount(&oxidation, 1.0, "SO3(g)"), 2.0);
    let hydration = MineralReaction::SulfurTrioxideHydration.build().unwrap();
    close(amount(&hydration, 2.0, "SO3(g)"), 2.0);
    close(amount(&hydration, 2.0, "H2O(l)"), 2.0);
    close(amount(&hydration, 2.0, "H2SO4(l)"), 2.0);
    assert!(
        hydration
            .terms()
            .iter()
            .all(|term| term.species != id("H+(aq)"))
    );
    let overall = MineralReaction::PyriteToHematiteAndSulfuricAcid
        .build()
        .unwrap();
    close(amount(&overall, 0.25, "FeS2(s)"), 1.0);
    close(amount(&overall, 0.25, "O2(g)"), 3.75);
    close(amount(&overall, 0.25, "H2O(l)"), 2.0);
    close(amount(&overall, 0.25, "Fe2O3(s)"), 0.5);
    close(amount(&overall, 0.25, "H2SO4(l)"), 2.0);
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

fn complete_inventory(
    reaction: &ConservedReaction,
    supplied: &[(&str, f64)],
) -> Vec<(SpeciesId, Amount)> {
    reaction
        .terms()
        .iter()
        .map(|term| {
            let amount = supplied
                .iter()
                .find_map(|(name, value)| (*name == term.species.as_str()).then_some(*value))
                .unwrap_or(0.0);
            (term.species.clone(), Amount::new(amount))
        })
        .collect()
}

fn inventory_amount(inventory: &[(SpeciesId, Amount)], species: &str) -> f64 {
    inventory
        .iter()
        .find(|(name, _)| name.as_str() == species)
        .unwrap()
        .1
        .value()
}

#[test]
fn bicarbonate_endpoints_retain_carbon_and_halve_theoretical_proton_demand() {
    for (acidic, bicarbonate, mineral, acid_moles, bicarbonate_moles) in [
        (
            MineralReaction::CalciteAcidNeutralization,
            MineralReaction::CalciteBicarbonateNeutralization,
            "CaCO3(s)",
            2.0,
            1.0,
        ),
        (
            MineralReaction::DolomiteAcidNeutralization,
            MineralReaction::DolomiteBicarbonateNeutralization,
            "CaMg(CO3)2(s)",
            4.0,
            2.0,
        ),
    ] {
        let acidic = acidic.build().unwrap();
        let bicarbonate = bicarbonate.build().unwrap();
        assert_eq!(amount(&acidic, 1.0, mineral), 1.0);
        assert_eq!(amount(&bicarbonate, 1.0, mineral), 1.0);
        assert_eq!(amount(&acidic, 1.0, "H+(aq)"), acid_moles);
        assert_eq!(amount(&bicarbonate, 1.0, "H+(aq)"), bicarbonate_moles);
        assert_eq!(amount(&bicarbonate, 1.0, "HCO3[-](aq)"), bicarbonate_moles);
        assert_ne!(
            acidic.conservation().identity(),
            bicarbonate.conservation().identity()
        );
        assert!(
            !bicarbonate
                .terms()
                .iter()
                .any(|term| term.species == id("CO2(aq)"))
        );
    }
}

#[test]
fn hydrated_lime_and_elsner_have_independent_known_answer_ratios() {
    let lime = MineralReaction::HydratedLimeNeutralization.build().unwrap();
    assert_eq!(amount(&lime, 10.0, "Ca(OH)2(s)"), 10.0);
    assert_eq!(amount(&lime, 10.0, "H+(aq)"), 20.0);
    assert_eq!(amount(&lime, 10.0, "H2O(l)"), 20.0);
    let gold = MineralReaction::GoldCyanidation.build().unwrap();
    assert_eq!(amount(&gold, 0.25, "Au(s)"), 1.0);
    assert_eq!(amount(&gold, 0.25, "CN[-](aq)"), 2.0);
    assert_eq!(amount(&gold, 0.25, "O2(g)"), 0.25);
    assert_eq!(amount(&gold, 0.25, "H2O(l)"), 0.5);
    assert_eq!(amount(&gold, 0.25, "Au(CN)2[-](aq)"), 1.0);
    assert_eq!(amount(&gold, 0.25, "OH[-](aq)"), 1.0);
}

#[test]
fn explicit_malachite_extent_consumes_protons_and_preserves_all_other_inventory() {
    let reaction = MineralReaction::MalachiteAcidLeach.build().unwrap();
    let initial = complete_inventory(
        &reaction,
        &[
            ("Cu2CO3(OH)2(s)", 10.0),
            ("H+(aq)", 12.0),
            ("H2O(l)", 100.0),
        ],
    );
    let original = initial.clone();
    let extent = reaction
        .limiting_extent(&[
            (id("Cu2CO3(OH)2(s)"), Amount::new(10.0)),
            (id("H+(aq)"), Amount::new(12.0)),
        ])
        .unwrap();
    assert_eq!(extent.value(), 3.0);
    let final_inventory = reaction.advance_inventory(&initial, extent).unwrap();
    assert_eq!(inventory_amount(&final_inventory, "Cu2CO3(OH)2(s)"), 7.0);
    assert_eq!(inventory_amount(&final_inventory, "H+(aq)"), 0.0);
    assert_eq!(inventory_amount(&final_inventory, "Cu[2+](aq)"), 6.0);
    assert_eq!(inventory_amount(&final_inventory, "CO2(aq)"), 3.0);
    assert_eq!(inventory_amount(&final_inventory, "H2O(l)"), 109.0);
    assert_eq!(initial, original);
    assert!(matches!(
        reaction.advance_inventory(&initial, Amount::new(3.001)),
        Err(ReactionBalanceError::ExtentExceedsInventory { .. })
    ));
    let mut reversed = initial;
    reversed.reverse();
    assert_eq!(
        reaction.advance_inventory(&reversed, extent).unwrap(),
        final_inventory
    );
    assert_eq!(
        reaction
            .advance_inventory(&reversed, Amount::new(0.0))
            .unwrap(),
        original
    );
}

#[test]
fn inventory_advance_preserves_partial_products_and_agrees_with_split_extents() {
    let reaction = MineralReaction::LimeHydration.build().unwrap();
    let initial = complete_inventory(
        &reaction,
        &[("CaO(s)", 10.0), ("H2O(l)", 12.0), ("Ca(OH)2(s)", 4.0)],
    );
    let intermediate = reaction
        .advance_inventory(&initial, Amount::new(2.0))
        .unwrap();
    let split = reaction
        .advance_inventory(&intermediate, Amount::new(3.0))
        .unwrap();
    let whole = reaction
        .advance_inventory(&initial, Amount::new(5.0))
        .unwrap();
    assert_eq!(split, whole);
    assert_eq!(inventory_amount(&whole, "Ca(OH)2(s)"), 9.0);
    assert_eq!(inventory_amount(&whole, "CaO(s)"), 5.0);
    assert_eq!(inventory_amount(&whole, "H2O(l)"), 7.0);
}

#[test]
fn complete_inventory_rejects_missing_duplicate_wrong_phase_and_invalid_amounts() {
    let reaction = MineralReaction::LimeHydration.build().unwrap();
    let good = complete_inventory(&reaction, &[("CaO(s)", 1.0), ("H2O(l)", 2.0)]);
    assert!(
        reaction
            .advance_inventory(&good[..2], Amount::new(0.5))
            .is_err()
    );
    let mut duplicate = good.clone();
    duplicate[1] = duplicate[0].clone();
    assert!(
        reaction
            .advance_inventory(&duplicate, Amount::new(0.5))
            .is_err()
    );
    let mut wrong_phase = good.clone();
    wrong_phase[0].0 = id("Ca(OH)2(aq)");
    assert!(
        reaction
            .advance_inventory(&wrong_phase, Amount::new(0.5))
            .is_err()
    );
    for invalid in [f64::NAN, f64::INFINITY, -1.0] {
        let mut bad = good.clone();
        bad[0].1 = Amount::new(invalid);
        assert!(reaction.advance_inventory(&bad, Amount::new(0.0)).is_err());
        assert!(
            reaction
                .advance_inventory(&good, Amount::new(invalid))
                .is_err()
        );
    }
}

#[test]
fn inventory_advance_refuses_vanishing_or_overflowing_changes_without_mutation() {
    let reaction = MineralReaction::LimeHydration.build().unwrap();
    let huge = complete_inventory(&reaction, &[("CaO(s)", 1e300), ("H2O(l)", 1e300)]);
    assert!(matches!(
        reaction.advance_inventory(&huge, Amount::new(1.0)),
        Err(ReactionBalanceError::Arithmetic(_))
    ));
    let overflow = complete_inventory(
        &reaction,
        &[
            ("CaO(s)", f64::MAX),
            ("H2O(l)", f64::MAX),
            ("Ca(OH)2(s)", f64::MAX),
        ],
    );
    let original = overflow.clone();
    assert!(matches!(
        reaction.advance_inventory(&overflow, Amount::new(f64::MAX)),
        Err(ReactionBalanceError::Arithmetic(_))
    ));
    assert_eq!(overflow, original);
    // Zero advancement is valid even when no positive extent is representable.
    let water = water([-2, -1, 2], [0, 0, 0]).unwrap();
    let tiny = complete_inventory(&water, &[("H2", f64::from_bits(1)), ("O2", 1.0)]);
    assert_eq!(
        water.advance_inventory(&tiny, Amount::new(0.0)).unwrap(),
        tiny
    );
}

#[test]
fn conservative_limiting_extent_does_not_hide_overdraw_in_rounded_products() {
    let reaction = MineralReaction::PyriteOxidationFerricHydroxide
        .build()
        .unwrap();
    // As exact binary fractions, fl(0.1 / 15) * 15 exceeds fl(0.1) by 2^-60,
    // but an ordinary multiplication rounds the product back onto stock.
    let nearest = 0.1_f64 / 15.0;
    assert_eq!(nearest * 15.0, 0.1);
    assert!(nearest.mul_add(15.0, -0.1) > 0.0);
    for available in [0.1, 1.0, 10.0, f64::MIN_POSITIVE, 1e200] {
        let reactants = [
            (id("FeS2(s)"), Amount::new(available)),
            (id("O2(g)"), Amount::new(available)),
            (id("H2O(l)"), Amount::new(available)),
        ];
        let limit = reaction.limiting_extent(&reactants).unwrap();
        assert!(limit.value().mul_add(15.0, -available) <= 0.0);
        assert!(limit.value().next_up().mul_add(15.0, -available) > 0.0);
        let initial = complete_inventory(
            &reaction,
            &[
                ("FeS2(s)", available),
                ("O2(g)", available),
                ("H2O(l)", available),
            ],
        );
        let advanced = reaction.advance_inventory(&initial, limit).unwrap();
        assert!(advanced.iter().all(|(_, value)| value.value() >= 0.0));
        assert!(matches!(
            reaction.advance_inventory(&initial, Amount::new(limit.value().next_up())),
            Err(ReactionBalanceError::ExtentExceedsInventory { .. })
        ));
    }
}
