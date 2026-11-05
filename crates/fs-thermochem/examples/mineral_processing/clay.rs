//! Source-bound copper/clay acid competition and subsequent aluminum hydrolysis.

use super::{NBS_REFERENCE, ReactionCase, load, reference_point};
use fs_matdb::{MatDbError, PropertyKey, QueryPoint, SelectionPolicy};
use fs_qty::{Amount, Dims, MolarMass, QuantitySpec};
use fs_thermochem::stoichiometry::{
    FormationEnthalpySource, minerals::MineralReaction, network::ReactionNetwork,
};
use std::{error::Error, path::Path};

pub(super) const USAGE: &str = "       mineral_processing --clay-acid-inventory PACK_DIR CUO_MOL KAOLINITE_MOL PROTON_MOL WATER_MOL CUO_FIRST_EXTENT_MOL KAOLINITE_EXTENT_MOL GIBBSITE_EXTENT_MOL CUO_LATER_EXTENT_MOL
Three explicit stages: shared-acid leaching, selected gibbsite hydrolysis,
then later CuO leaching. Products initially zero. Extents are input decisions,
not rates or equilibrium predictions. Gibbsite extent is per Al2O3.3H2O unit.
See data/matdb/seed-v1/MINING_CLAY_ACID_DEMAND.md for the eight required packs.";

pub(super) fn case(name: &str) -> Option<ReactionCase> {
    Some(match name {
        "kaolinite-acid-leach" => ReactionCase {
            reaction: MineralReaction::KaoliniteAcidLeach,
            reference: NBS_REFERENCE,
            packs: &[
                "kaolinite-solid-nbs1989-thermochemistry",
                "hydrogen-ion-aqueous-nbs1989",
                "aluminum-iii-ion-aqueous-nbs1989",
                "silicic-acid-component-aqueous-nbs1989",
                "water-liquid-nbs1989-thermochemistry",
            ],
            species: &[
                "Al2Si2O5(OH)4(s-kaolinite)",
                "H+(aq)",
                "Al[3+](aq)",
                "H4SiO4(aq)",
                "H2O(l)",
            ],
            phase_axes: &[
                "source-phase-solid",
                "source-phase-aqueous",
                "source-phase-aqueous",
                "source-phase-aqueous",
                "source-phase-liquid",
            ],
            equation: "Al2Si2O5(OH)4(kaolinite) + 6 H+(aq) -> 2 Al[3+](aq) + 2 H4SiO4(aq) + H2O(l), retaining the NBS formal silicic-acid component",
        },
        "aluminum-hydrolysis-gibbsite" => ReactionCase {
            reaction: MineralReaction::AluminumHydrolysisGibbsite,
            reference: NBS_REFERENCE,
            packs: &[
                "aluminum-iii-ion-aqueous-nbs1989",
                "water-liquid-nbs1989-thermochemistry",
                "gibbsite-solid-nbs1989-thermochemistry",
                "hydrogen-ion-aqueous-nbs1989",
            ],
            species: &["Al[3+](aq)", "H2O(l)", "Al2O3.3H2O(s-gibbsite)", "H+(aq)"],
            phase_axes: &[
                "source-phase-aqueous",
                "source-phase-liquid",
                "source-phase-solid",
                "source-phase-aqueous",
            ],
            equation: "2 Al[3+](aq) + 6 H2O(l) -> Al2O3.3H2O(gibbsite) + 6 H+(aq), one doubled NBS gibbsite unit per extent",
        },
        "cupric-oxide-acid-leach-nbs1989" => ReactionCase {
            reaction: MineralReaction::CupricOxideAcidLeach,
            reference: NBS_REFERENCE,
            packs: &[
                "copper-ii-oxide-solid-nbs1989-thermochemistry",
                "hydrogen-ion-aqueous-nbs1989",
                "copper-ii-ion-aqueous-nbs1989",
                "water-liquid-nbs1989-thermochemistry",
            ],
            species: &["CuO(s)", "H+(aq)", "Cu[2+](aq)", "H2O(l)"],
            phase_axes: &[
                "source-phase-solid",
                "source-phase-aqueous",
                "source-phase-aqueous",
                "source-phase-liquid",
            ],
            equation: "CuO(cr) + 2 H+(aq) -> Cu[2+](aq) + H2O(l), matched NBS source family",
        },
        _ => return None,
    })
}

pub(super) fn qualify_point(mut at: QueryPoint, pack: &str) -> Result<QueryPoint, MatDbError> {
    let flags: &[(&str, f64)] = match pack {
        "kaolinite-solid-nbs1989-thermochemistry" => &[("source-phase-kaolinite", 1.0)],
        "gibbsite-solid-nbs1989-thermochemistry" => &[
            ("source-phase-gibbsite", 1.0),
            ("source-formula-aluminum-atoms", 2.0),
        ],
        "silicic-acid-component-aqueous-nbs1989" => {
            &[("source-silicic-component-h2sio3-plus-water", 1.0)]
        }
        "copper-ii-oxide-solid-nbs1989-thermochemistry" => {
            &[("source-phase-copper-ii-oxide-crystalline", 1.0)]
        }
        _ => &[],
    };
    for &(axis, value) in flags {
        at = at.with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)?;
    }
    Ok(at)
}

pub(super) fn run(args: &[String]) -> Result<(), Box<dyn Error>> {
    if args.len() != 10 {
        return Err(USAGE.into());
    }
    let values = args[2..]
        .iter()
        .map(|v| v.parse::<f64>())
        .collect::<Result<Vec<_>, _>>()?;
    let [
        cuo,
        clay,
        acid,
        water,
        first_cu,
        dissolve,
        precipitate,
        later_cu,
    ] = values[..]
    else {
        return Err(USAGE.into());
    };
    let reactions = [
        MineralReaction::CupricOxideAcidLeach,
        MineralReaction::KaoliniteAcidLeach,
        MineralReaction::AluminumHydrolysisGibbsite,
    ]
    .into_iter()
    .map(|r| r.build())
    .collect::<Result<Vec<_>, _>>()?;
    let ids: Vec<_> = reactions.iter().map(|r| r.id().clone()).collect();
    let network = ReactionNetwork::new(reactions)?;
    let initial: Vec<_> = network
        .species()
        .iter()
        .map(|id| {
            let n = match id.as_str() {
                "CuO(s)" => cuo,
                "Al2Si2O5(OH)4(s-kaolinite)" => clay,
                "H+(aq)" => acid,
                "H2O(l)" => water,
                _ => 0.0,
            };
            (id.clone(), Amount::new(n))
        })
        .collect();
    // Reuse the exact source associations and state selectors of the single-
    // reaction modes. Each shared species has one source in every stage.
    let cases = [
        "cupric-oxide-acid-leach-nbs1989",
        "kaolinite-acid-leach",
        "aluminum-hydrolysis-gibbsite",
    ]
    .map(|name| case(name).ok_or("missing clay source case"))
    .into_iter()
    .collect::<Result<Vec<_>, _>>()?;
    let sources = network
        .species()
        .iter()
        .map(|id| {
            let (selected, index) = cases
                .iter()
                .find_map(|c| {
                    c.species
                        .iter()
                        .position(|s| *s == id.as_str())
                        .map(|i| (c, i))
                })
                .ok_or("missing species source association")?;
            let slug = selected.packs[index];
            let pack = load(&Path::new(&args[1]).join(format!("{slug}.fsmatpk")), slug)?;
            let point = reference_point(NBS_REFERENCE, selected.phase_axes[index], slug)?;
            let mass = pack.claims().query_typed(
                &PropertyKey::new("molar-mass", MolarMass::DIMS),
                &point,
                SelectionPolicy::SingleClaimOnly,
            )?;
            pack.claims().verify_receipt(&mass.receipt)?;
            Ok((id.clone(), pack, point, mass.evidence.value.value))
        })
        .collect::<Result<Vec<_>, Box<dyn Error>>>()?;
    let heat: Vec<_> = sources
        .iter()
        .map(|(id, pack, at, _)| FormationEnthalpySource {
            species: id,
            pack,
            at,
            reference_axis: NBS_REFERENCE,
        })
        .collect();
    let mut current = initial.clone();
    let mut stages = Vec::new();
    for amounts in [
        [first_cu, dissolve, 0.0],
        [0.0, 0.0, precipitate],
        [later_cu, 0.0, 0.0],
    ] {
        let requested: Vec<_> = ids.iter().cloned().zip(amounts.map(Amount::new)).collect();
        let stage = network.advance_with_standard_enthalpy(&current, &requested, &heat)?;
        current = stage.inventory.final_inventory.clone();
        stages.push(stage);
    }
    let mass = |inventory: &[(fs_thermochem::SpeciesId, Amount)]| {
        inventory
            .iter()
            .zip(&sources)
            .map(|((_, n), (_, _, _, m))| n.value() * m)
            .sum::<f64>()
    };
    let mass_residual = mass(&current) - mass(&initial);
    let enthalpy: f64 = stages.iter().map(|s| s.enthalpy.value()).sum();
    if !mass_residual.is_finite() || !enthalpy.is_finite() {
        return Err("unrepresentable source-mass or standard-enthalpy total".into());
    }
    // Publish only after all three finite stages and source queries succeed.
    println!("Prescribed copper/clay acid accounting; NBS 298.15 K, 100000 Pa.");
    println!("stage,reaction,extent_mol,standard_delta_h_J");
    for (name, stage) in [
        "shared-acid-leach",
        "gibbsite-hydrolysis",
        "later-copper-leach",
    ]
    .into_iter()
    .zip(&stages)
    {
        for reaction in &stage.reactions {
            println!(
                "{name},{},{},{}",
                reaction.reaction.as_str(),
                reaction.extent.value(),
                reaction.enthalpy.value()
            );
        }
    }
    println!(
        "species,initial_mol,after_shared_leach_mol,after_hydrolysis_mol,final_mol,final_species_mass_kg"
    );
    for (i, (id, _, _, m)) in sources.iter().enumerate() {
        println!(
            "{},{},{},{},{},{}",
            id.as_str(),
            initial[i].1.value(),
            stages[0].inventory.final_inventory[i].1.value(),
            stages[1].inventory.final_inventory[i].1.value(),
            current[i].1.value(),
            current[i].1.value() * m
        );
    }
    println!("Total standard chemical enthalpy change [J]: {enthalpy}");
    println!("Source-molar-mass residual, final minus initial [kg]: {mass_residual}");
    for (id, pack, _, _) in &sources {
        println!(
            "source,{},{},{:?}",
            id.as_str(),
            pack.pack_id(),
            pack.content_hash()
        );
    }
    println!(
        "Atoms and charge close exactly; independently rounded source molar masses are retained. Proton inventory is not a commercial acid dose. Gibbsite uses Al2O3.3H2O, and H4SiO4 retains its formal water component. Extents and phase selection are declared: no clay dissolution rate, pH, nonideal speciation, automatic precipitation, recovery or operating heat duty is inferred."
    );
    Ok(())
}
