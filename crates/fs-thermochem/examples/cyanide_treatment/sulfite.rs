//! Staged source-pack consumer for the free-cyanide SO2/air endpoint.

use super::{heat_sources, load, receipts};
use fs_qty::Amount;
use fs_thermochem::stoichiometry::{minerals::MineralReaction, network::ReactionNetwork};
use std::{error::Error, path::Path};

pub(super) fn run(directory: &Path, values: &[f64]) -> Result<(), Box<dyn Error>> {
    let [cn, reagent, water, oxygen, dissolve, productive, side] = *values else {
        return Err("sulfite requires four initial amounts and three reaction extents".into());
    };
    let reactions = [
        MineralReaction::SodiumMetabisulfiteDissolution,
        MineralReaction::CyanideBisulfiteOxidation,
        MineralReaction::BisulfiteOxygenOxidation,
    ]
    .into_iter()
    .map(|r| r.build())
    .collect::<Result<Vec<_>, _>>()?;
    let requested = |amounts: [f64; 3]| {
        reactions
            .iter()
            .zip(amounts)
            .map(|(r, n)| (r.id().clone(), Amount::new(n)))
            .collect::<Vec<_>>()
    };
    let preparation = requested([dissolve, 0.0, 0.0]);
    let oxidation = requested([0.0, productive, side]);
    let network = ReactionNetwork::new(reactions)?;
    let stocks = [
        ("CN[-](aq)", cn),
        ("Na2S2O5(s)", reagent),
        ("H2O(l)", water),
        ("O2(aq)", oxygen),
    ];
    let initial: Vec<_> = network
        .species()
        .iter()
        .map(|id| {
            let n = stocks
                .iter()
                .find(|(name, _)| *name == id.as_str())
                .map_or(0.0, |(_, n)| *n);
            (id.clone(), Amount::new(n))
        })
        .collect();
    let sources = network
        .species()
        .iter()
        .map(|s| load(directory, s))
        .collect::<Result<Vec<_>, _>>()?;
    let heat = heat_sources(&sources);
    let prepared = network.advance_with_standard_enthalpy(&initial, &preparation, &heat)?;
    let treated = network.advance_with_standard_enthalpy(
        &prepared.inventory.final_inventory,
        &oxidation,
        &heat,
    )?;

    // Resolve all stages and diagnostics before publishing an accepted result.
    let mass = |inventory: &[(fs_thermochem::SpeciesId, Amount)]| {
        inventory
            .iter()
            .zip(&sources)
            .map(|((_, n), s)| n.value() * s.molar_mass)
            .sum::<f64>()
    };
    let residual = mass(&treated.inventory.final_inventory) - mass(&initial);
    let enthalpy = prepared.enthalpy.value() + treated.enthalpy.value();
    if !residual.is_finite() || !enthalpy.is_finite() {
        return Err("unrepresentable source-mass or enthalpy total".into());
    }
    println!("Prescribed finite sulfite treatment; NBS 298.15 K, 100000 Pa.");
    println!("stage,reaction,extent_mol,standard_delta_h_J");
    for (stage_name, stage) in [
        ("reagent-preparation", &prepared),
        ("productive-and-side-oxidation", &treated),
    ] {
        for reaction in &stage.reactions {
            println!(
                "{stage_name},{},{},{}",
                reaction.reaction.as_str(),
                reaction.extent.value(),
                reaction.enthalpy.value()
            );
        }
    }
    println!("species,initial_mol,after_preparation_mol,after_treatment_mol,final_species_mass_kg");
    for (i, s) in sources.iter().enumerate() {
        let final_amount = treated.inventory.final_inventory[i].1.value();
        println!(
            "{},{},{},{},{}",
            s.species.as_str(),
            initial[i].1.value(),
            prepared.inventory.final_inventory[i].1.value(),
            final_amount,
            final_amount * s.molar_mass
        );
    }
    println!("Total standard chemical enthalpy change [J]: {enthalpy}");
    println!("Source-mass residual, final minus initial [kg]: {residual}");
    receipts(&sources);
    println!(
        "Free cyanide and cyanate retain nitrogen; sulfate, sodium and proton equivalents remain explicit. No catalyst kinetics, WAD-CN conversion, automatic pH, gas transfer, plant dose or operating heat duty is inferred. Na2S2O5(s) has no supplied Gibbs/entropy data."
    );
    Ok(())
}
