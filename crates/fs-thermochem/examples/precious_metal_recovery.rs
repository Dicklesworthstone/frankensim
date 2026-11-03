//! Finite native-silver cyanidation and shared-zinc gold/silver cementation.
//!
//! Compile the named seed packs with `xtask matdb-pack` into PACK_DIR; each
//! filename is PACK_ID.fsmatpk. Run with --help for the two explicit stage
//! inputs. The existing reaction-network owner admits all gross reagent
//! demands before returning any material or standard-enthalpy result.
//!
//! Extents are caller decisions, not a rate, equilibrium, recovery allocation,
//! or plant calibration. The source state is 298.15 K and 100000 Pa. Aqueous
//! stocks are moles of the named species, not solution volume or reagent mass.

use std::{error::Error, path::Path};

use fs_matdb::{NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy};
use fs_qty::{
    Amount, Dims, MolarMass, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{
    SpeciesId,
    stoichiometry::{
        FormationEnthalpySource, minerals::MineralReaction, network::ReactionNetwork,
        thermodynamics::FormationGibbsSource,
    },
};

const REFERENCE: &str = "source-elemental-reference-nbs1989";
const USAGE: &str = "usage:
  precious_metal_recovery silver-leach PACK_DIR AG_MOL CN_MOL DISSOLVED_O2_MOL WATER_MOL EXTENT_MOL
  precious_metal_recovery cementation PACK_DIR AU_COMPLEX_MOL AG_COMPLEX_MOL ZN_MOL AU_EXTENT_MOL AG_EXTENT_MOL

Each extent is moles of the complete integer equation as written:
  4 Ag + 8 CN- + O2(aq) + 2 H2O -> 4 Ag(CN)2- + 4 OH-
  2 Au(CN)2- + Zn -> 2 Au + Zn(CN)4(2-)
  2 Ag(CN)2- + Zn -> 2 Ag + Zn(CN)4(2-)

All initial products are zero. Both cementation pathways share one zinc stock.
No automatic recovery allocation, side-reaction demand, or kinetics is supplied.
PACK_DIR must contain the phase-specific NBS .fsmatpk files documented in
data/matdb/seed-v1/MINING_PRECIOUS_METALS.md.";

struct Source {
    species: SpeciesId,
    pack: NormalizedPack,
    at: QueryPoint,
    molar_mass: f64,
}

fn point(phase: &str, ion: bool) -> Result<QueryPoint, Box<dyn Error>> {
    let semantic = |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
    let mut at = QueryPoint::new()
        .with_quantity(
            "temperature",
            semantic(QuantityKind::AbsoluteTemperature),
            298.15,
        )?
        .with_quantity("pressure", semantic(QuantityKind::Pressure), 100_000.0)?;
    let mut flags = vec![
        (REFERENCE, 1.0),
        ("source-standard-state", 1.0),
        (phase, 1.0),
        ("engineering-reference-model", 1.0),
        ("source-specimen-qualified", 0.0),
    ];
    if phase == "source-phase-aqueous" {
        flags.push(("source-aqueous-unit-activity-infinite-dilution", 1.0));
        if ion {
            flags.extend([
                ("source-individual-ion", 1.0),
                ("source-ionic-reference-hydrogen-zero", 1.0),
            ]);
        }
    } else {
        flags.push(("source-pure-species", 1.0));
    }
    for (axis, value) in flags {
        at = at.with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)?;
    }
    Ok(at)
}

fn load(directory: &Path, species: &SpeciesId) -> Result<Source, Box<dyn Error>> {
    let (slug, phase, ion) = match species.as_str() {
        "Ag(s)" => (
            "silver-solid-nbs1989-thermochemistry",
            "source-phase-solid",
            false,
        ),
        "Au(s)" => (
            "gold-solid-nbs1989-thermochemistry",
            "source-phase-solid",
            false,
        ),
        "Zn(s)" => (
            "zinc-solid-nbs1989-thermochemistry",
            "source-phase-solid",
            false,
        ),
        "Ag(CN)2[-](aq)" => (
            "dicyanoargentate-ion-aqueous-nbs1989",
            "source-phase-aqueous",
            true,
        ),
        "Au(CN)2[-](aq)" => (
            "dicyanoaurate-ion-aqueous-nbs1989",
            "source-phase-aqueous",
            true,
        ),
        "Zn(CN)4[2-](aq)" => (
            "tetracyanozincate-ion-aqueous-nbs1989",
            "source-phase-aqueous",
            true,
        ),
        "CN[-](aq)" => ("cyanide-ion-aqueous-nbs1989", "source-phase-aqueous", true),
        "OH[-](aq)" => (
            "hydroxide-ion-aqueous-nbs1989",
            "source-phase-aqueous",
            true,
        ),
        "O2(aq)" => (
            "oxygen-aqueous-nbs1989-thermochemistry",
            "source-phase-aqueous",
            false,
        ),
        "H2O(l)" => (
            "water-liquid-nbs1989-thermochemistry",
            "source-phase-liquid",
            false,
        ),
        _ => return Err(format!("no explicit NBS source association for {species:?}").into()),
    };
    let path = directory.join(format!("{slug}.fsmatpk"));
    let bytes = std::fs::read(&path).map_err(|error| format!("{}: {error}", path.display()))?;
    let pack = NormalizedPack::from_bytes(&bytes)?;
    if pack.pack_id() != slug {
        return Err(format!("{} must contain pack {slug}", path.display()).into());
    }
    let at = point(phase, ion)?;
    let mass = pack.claims().query_typed(
        &PropertyKey::new("molar-mass", MolarMass::DIMS),
        &at,
        SelectionPolicy::SingleClaimOnly,
    )?;
    pack.claims().verify_receipt(&mass.receipt)?;
    let molar_mass = mass.evidence.value.value;
    if !molar_mass.is_finite() || molar_mass <= 0.0 {
        return Err("source molar mass must be finite and positive".into());
    }
    Ok(Source {
        species: species.clone(),
        pack,
        at,
        molar_mass,
    })
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() == 1 && matches!(args[0].as_str(), "--help" | "-h") {
        println!("{USAGE}");
        return Ok(());
    }
    if args.len() != 7 {
        return Err(USAGE.into());
    }
    let values = args[2..]
        .iter()
        .map(|text| text.parse::<f64>())
        .collect::<Result<Vec<_>, _>>()?;
    let (recipes, stocks, extents) = match args[0].as_str() {
        "silver-leach" => (
            vec![MineralReaction::SilverCyanidationDissolvedOxygen],
            vec![
                ("Ag(s)", values[0]),
                ("CN[-](aq)", values[1]),
                ("O2(aq)", values[2]),
                ("H2O(l)", values[3]),
            ],
            vec![values[4]],
        ),
        "cementation" => (
            vec![
                MineralReaction::GoldCyanideCementationWithZinc,
                MineralReaction::SilverCyanideCementationWithZinc,
            ],
            vec![
                ("Au(CN)2[-](aq)", values[0]),
                ("Ag(CN)2[-](aq)", values[1]),
                ("Zn(s)", values[2]),
            ],
            vec![values[3], values[4]],
        ),
        _ => return Err(USAGE.into()),
    };
    let reactions = recipes
        .into_iter()
        .map(|recipe| recipe.build())
        .collect::<Result<Vec<_>, _>>()?;
    let requested: Vec<_> = reactions
        .iter()
        .zip(extents)
        .map(|(reaction, extent)| (reaction.id().clone(), Amount::new(extent)))
        .collect();
    let network = ReactionNetwork::new(reactions)?;
    let initial: Vec<_> = network
        .species()
        .iter()
        .map(|species| {
            let amount = stocks
                .iter()
                .find(|(name, _)| *name == species.as_str())
                .map_or(0.0, |(_, amount)| *amount);
            (species.clone(), Amount::new(amount))
        })
        .collect();
    let sources = network
        .species()
        .iter()
        .map(|species| load(Path::new(&args[1]), species))
        .collect::<Result<Vec<_>, _>>()?;
    let heat_sources: Vec<_> = sources
        .iter()
        .map(|source| FormationEnthalpySource {
            species: &source.species,
            pack: &source.pack,
            at: &source.at,
            reference_axis: REFERENCE,
        })
        .collect();
    let result = network.advance_with_standard_enthalpy(&initial, &requested, &heat_sources)?;

    // Resolve every requested diagnostic before printing an accepted result.
    let mut diagnostics = Vec::new();
    for reaction in network.reactions() {
        let active: Vec<_> = reaction
            .terms()
            .iter()
            .map(|term| {
                sources
                    .iter()
                    .find(|source| source.species == term.species)
                    .ok_or("missing active source")
            })
            .collect::<Result<Vec<_>, _>>()?;
        let h: Vec<_> = active
            .iter()
            .map(|s| FormationEnthalpySource {
                species: &s.species,
                pack: &s.pack,
                at: &s.at,
                reference_axis: REFERENCE,
            })
            .collect();
        let g: Vec<_> = active
            .iter()
            .map(|s| FormationGibbsSource {
                species: &s.species,
                pack: &s.pack,
                at: &s.at,
                reference_axis: REFERENCE,
            })
            .collect();
        let caloric = reaction.standard_equilibrium(&h)?;
        let direct = reaction.standard_gibbs_change(&g)?;
        diagnostics.push((
            reaction.id(),
            caloric.molar_gibbs().value(),
            direct.molar().value(),
        ));
    }
    println!(
        "Finite prescribed-extent {} stage; NBS reference at 298.15 K and 100000 Pa.",
        args[0]
    );
    println!(
        "reaction,extent_mol,standard_delta_h_J,standard_delta_g_from_h_s_J_per_mol,direct_standard_delta_g_J_per_mol"
    );
    for reaction in &result.reactions {
        let (_, caloric, direct) = diagnostics
            .iter()
            .find(|(id, _, _)| **id == reaction.reaction)
            .ok_or("missing reaction diagnostic")?;
        println!(
            "{},{},{},{},{}",
            reaction.reaction.as_str(),
            reaction.extent.value(),
            reaction.enthalpy.value(),
            caloric,
            direct
        );
    }
    println!(
        "Total standard chemical enthalpy change [J]: {}",
        result.enthalpy.value()
    );
    println!("species,initial_mol,consumed_mol,formed_mol,final_mol,final_species_mass_kg");
    for (index, source) in sources.iter().enumerate() {
        // The original network provides the same complete canonical species axis.
        let final_mol = result.inventory.final_inventory[index].1.value();
        println!(
            "{},{},{},{},{},{}",
            source.species.as_str(),
            initial[index].1.value(),
            result.inventory.consumed[index].1.value(),
            result.inventory.formed[index].1.value(),
            final_mol,
            final_mol * source.molar_mass
        );
    }
    println!("source_species,pack_id,content_hash");
    for source in sources {
        println!(
            "{},{},{:?}",
            source.species.as_str(),
            source.pack.pack_id(),
            source.pack.content_hash()
        );
    }
    println!(
        "Standard-state reaction accounting only. No kinetic recovery allocation, side reactions, solution activity model, deaeration, filtration, heat of mixing, or operating-temperature duty is inferred."
    );
    Ok(())
}
