//! Prescribed gold/silver-sulfide leaching and sulfur cyanide demand.
//!
//! The first stage shares finite CN-, dissolved oxygen and water across gold,
//! silver-sulfide and elemental-sulfur endpoints. The second stage transfers
//! sulfur from the resulting thiosulfate to the remaining cyanide. This explicit
//! ordering is an accounting input, not a kinetic or universal mechanism claim.
//! Compile the twelve NBS packs listed in MINING_SULFUR_CYANIDE.md into PACK_DIR.

use std::{error::Error, path::Path};

use fs_matdb::{NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy};
use fs_qty::{
    Amount, Dims, MolarMass, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{
    ReactionId, SpeciesId,
    stoichiometry::{
        FormationEnthalpySource, minerals::MineralReaction, network::ReactionNetwork,
        thermodynamics::FormationGibbsSource,
    },
};

const REFERENCE: &str = "source-elemental-reference-nbs1989";
const AG2S: &str = "Ag2S(s-nbs-alpha-orthorhombic)";
const USAGE: &str = "usage:
  sulfur_cyanide_demand PACK_DIR AU_MOL AG2S_MOL S_MOL CN_MOL O2_AQ_MOL WATER_MOL AU_EXTENT_MOL AG2S_EXTENT_MOL S_EXTENT_MOL THIOSULFATE_EXTENT_MOL

Initial products are explicitly zero. Amounts are moles of the named species.
Each extent is moles of the complete equation as written:
  4 Au + 8 CN- + O2(aq) + 2 H2O -> 4 Au(CN)2- + 4 OH-
  2 Ag2S + 8 CN- + 2 O2(aq) + H2O -> 4 Ag(CN)2- + S2O3(2-) + 2 OH-
  S(s) + CN- -> SCN-
  S2O3(2-) + CN- -> SCN- + SO3(2-)

The first three extents share the original finite stocks. The fourth uses the
first stage's final inventory. S is per mole sulfur atoms, not per S8 molecule.
No recovery allocation, sulfur-product selection, pH or rate is inferred.
PACK_DIR contains the phase-specific NBS .fsmatpk files listed in
data/matdb/seed-v1/MINING_SULFUR_CYANIDE.md.";

struct Source {
    species: SpeciesId,
    pack: NormalizedPack,
    at: QueryPoint,
    molar_mass: f64,
}

fn load(directory: &Path, species: &SpeciesId) -> Result<Source, Box<dyn Error>> {
    let (slug, phase, ion) = match species.as_str() {
        AG2S => ("silver-sulfide-alpha-solid-nbs1989", "solid", false),
        "S2O3[2-](aq)" => ("thiosulfate-ion-aqueous-nbs1989", "aqueous", true),
        "SCN[-](aq)" => ("thiocyanate-ion-aqueous-nbs1989", "aqueous", true),
        "S(s)" => ("sulfur-orthorhombic-solid-nbs1989", "solid", false),
        "SO3[2-](aq)" => ("sulfite-ion-aqueous-nbs1989", "aqueous", true),
        "Au(s)" => ("gold-solid-nbs1989-thermochemistry", "solid", false),
        "Au(CN)2[-](aq)" => ("dicyanoaurate-ion-aqueous-nbs1989", "aqueous", true),
        "Ag(CN)2[-](aq)" => ("dicyanoargentate-ion-aqueous-nbs1989", "aqueous", true),
        "CN[-](aq)" => ("cyanide-ion-aqueous-nbs1989", "aqueous", true),
        "OH[-](aq)" => ("hydroxide-ion-aqueous-nbs1989", "aqueous", true),
        "O2(aq)" => ("oxygen-aqueous-nbs1989-thermochemistry", "aqueous", false),
        "H2O(l)" => ("water-liquid-nbs1989-thermochemistry", "liquid", false),
        _ => return Err(format!("no explicit NBS source association for {species:?}").into()),
    };
    let path = directory.join(format!("{slug}.fsmatpk"));
    let bytes = std::fs::read(&path).map_err(|error| format!("{}: {error}", path.display()))?;
    let pack = NormalizedPack::from_bytes(&bytes)?;
    if pack.pack_id() != slug {
        return Err(format!("{} must contain pack {slug}", path.display()).into());
    }
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
        ("engineering-reference-model", 1.0),
        ("source-specimen-qualified", 0.0),
        (
            match phase {
                "solid" => "source-phase-solid",
                "liquid" => "source-phase-liquid",
                _ => "source-phase-aqueous",
            },
            1.0,
        ),
        (
            if phase == "aqueous" {
                "source-aqueous-unit-activity-infinite-dilution"
            } else {
                "source-pure-species"
            },
            1.0,
        ),
    ];
    if ion {
        flags.extend([
            ("source-individual-ion", 1.0),
            ("source-ionic-reference-hydrogen-zero", 1.0),
        ]);
    }
    match species.as_str() {
        AG2S => flags.push(("source-phase-silver-sulfide-nbs-alpha-orthorhombic", 1.0)),
        "S(s)" => flags.push(("source-phase-sulfur-orthorhombic", 1.0)),
        _ => (),
    }
    for (axis, value) in flags {
        at = at.with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)?;
    }
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

fn requests(
    gold: f64,
    silver: f64,
    sulfur: f64,
    thiosulfate: f64,
) -> Result<Vec<(ReactionId, Amount)>, Box<dyn Error>> {
    [
        ("gold-cyanidation-elsner-dissolved-oxygen", gold),
        ("silver-sulfide-cyanidation-thiosulfate", silver),
        ("elemental-sulfur-cyanide-consumption", sulfur),
        ("thiosulfate-cyanide-sulfur-transfer", thiosulfate),
    ]
    .into_iter()
    .map(|(id, amount)| Ok((ReactionId::new(id)?, Amount::new(amount))))
    .collect()
}

fn amount(stock: &[(SpeciesId, Amount)], name: &str) -> Result<f64, Box<dyn Error>> {
    Ok(stock
        .iter()
        .find(|(id, _)| id.as_str() == name)
        .ok_or("missing admitted species")?
        .1
        .value())
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() == 1 && matches!(args[0].as_str(), "--help" | "-h") {
        println!("{USAGE}");
        return Ok(());
    }
    if args.len() != 11 {
        return Err(USAGE.into());
    }
    let values = args[1..]
        .iter()
        .map(|text| text.parse::<f64>())
        .collect::<Result<Vec<_>, _>>()?;
    if values.iter().any(|v| !v.is_finite() || *v < 0.0) {
        return Err("all species amounts and extents must be finite and nonnegative".into());
    }
    let reactions = [
        MineralReaction::GoldCyanidationDissolvedOxygen,
        MineralReaction::SilverSulfideCyanidationThiosulfate,
        MineralReaction::ElementalSulfurCyanideConsumption,
        MineralReaction::ThiosulfateCyanideSulfurTransfer,
    ]
    .into_iter()
    .map(|r| r.build())
    .collect::<Result<Vec<_>, _>>()?;
    let network = ReactionNetwork::new(reactions)?;
    let initial: Vec<_> = network
        .species()
        .iter()
        .map(|id| {
            let amount = match id.as_str() {
                "Au(s)" => values[0],
                AG2S => values[1],
                "S(s)" => values[2],
                "CN[-](aq)" => values[3],
                "O2(aq)" => values[4],
                "H2O(l)" => values[5],
                _ => 0.0, // The CLI explicitly declares every initial product zero.
            };
            (id.clone(), Amount::new(amount))
        })
        .collect();
    let sources = network
        .species()
        .iter()
        .map(|id| load(Path::new(&args[0]), id))
        .collect::<Result<Vec<_>, _>>()?;
    let heat: Vec<_> = sources
        .iter()
        .map(|s| FormationEnthalpySource {
            species: &s.species,
            pack: &s.pack,
            at: &s.at,
            reference_axis: REFERENCE,
        })
        .collect();
    let first = network.advance_with_standard_enthalpy(
        &initial,
        &requests(values[6], values[7], values[8], 0.0)?,
        &heat,
    )?;
    let second = network.advance_with_standard_enthalpy(
        &first.inventory.final_inventory,
        &requests(0.0, 0.0, 0.0, values[9])?,
        &heat,
    )?;
    // Admit both stages and every source diagnostic before printing a result.
    let mut diagnostics = Vec::new();
    for reaction in network.reactions() {
        let active = reaction
            .terms()
            .iter()
            .map(|term| {
                sources
                    .iter()
                    .find(|s| s.species == term.species)
                    .ok_or("missing reaction source")
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
        diagnostics.push((
            reaction.id().as_str().to_owned(),
            reaction.standard_equilibrium(&h)?.molar_gibbs().value(),
            reaction.standard_gibbs_change(&g)?.molar().value(),
        ));
    }
    let total_enthalpy = first.enthalpy.value() + second.enthalpy.value();
    if !total_enthalpy.is_finite() {
        return Err("total standard chemical enthalpy exceeds numerical range".into());
    }
    println!("Explicit sulfur cyanide-demand stages; NBS 298.15 K, 100000 Pa.");
    println!(
        "stage,reaction,extent_mol,standard_delta_h_J,delta_g_from_h_s_J_per_mol,direct_standard_delta_g_J_per_mol"
    );
    for (label, result) in [
        ("leach-and-sulfur", &first),
        ("thiosulfate-transfer", &second),
    ] {
        for reaction in &result.reactions {
            let (_, hs, direct) = diagnostics
                .iter()
                .find(|(id, _, _)| id == reaction.reaction.as_str())
                .ok_or("missing diagnostic")?;
            println!(
                "{label},{},{},{},{hs},{direct}",
                reaction.reaction.as_str(),
                reaction.extent.value(),
                reaction.enthalpy.value()
            );
        }
    }
    println!("stage,species,initial_mol,consumed_mol,formed_mol,final_mol,final_species_mass_kg");
    for (label, start, result) in [
        ("leach-and-sulfur", &initial, &first),
        (
            "thiosulfate-transfer",
            &first.inventory.final_inventory,
            &second,
        ),
    ] {
        for (i, source) in sources.iter().enumerate() {
            let n = result.inventory.final_inventory[i].1.value();
            println!(
                "{label},{},{},{},{},{n},{}",
                source.species.as_str(),
                start[i].1.value(),
                result.inventory.consumed[i].1.value(),
                result.inventory.formed[i].1.value(),
                n * source.molar_mass
            );
        }
    }
    println!(
        "stage,total_source_mass_kg,Au_atoms_mol,Ag_atoms_mol,S_atoms_mol,CN_units_mol,free_CN_mol"
    );
    for (label, stock) in [
        ("initial", &initial),
        ("after-leach", &first.inventory.final_inventory),
        ("after-thiosulfate", &second.inventory.final_inventory),
    ] {
        let n = |name| amount(stock, name);
        let mass: f64 = sources
            .iter()
            .zip(stock)
            .map(|(s, (_, n))| s.molar_mass * n.value())
            .sum();
        let gold = n("Au(s)")? + n("Au(CN)2[-](aq)")?;
        let silver = 2.0 * n(AG2S)? + n("Ag(CN)2[-](aq)")?;
        let sulfur =
            n(AG2S)? + n("S(s)")? + 2.0 * n("S2O3[2-](aq)")? + n("SCN[-](aq)")? + n("SO3[2-](aq)")?;
        let cn_units = n("CN[-](aq)")?
            + 2.0 * n("Au(CN)2[-](aq)")?
            + 2.0 * n("Ag(CN)2[-](aq)")?
            + n("SCN[-](aq)")?;
        println!(
            "{label},{mass},{gold},{silver},{sulfur},{cn_units},{}",
            n("CN[-](aq)")?
        );
    }
    println!("Total standard chemical enthalpy change [J]: {total_enthalpy}");
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
        "Prescribed extents and source-reference chemical energies only. No rates, pH, nonideal activities, heat of mixing, operating-temperature duty, ore accessibility or recovery forecast is inferred. Net charge is conserved for the modeled reactive species; spectator counterions are outside this inventory."
    );
    Ok(())
}
