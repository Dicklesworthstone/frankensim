//! Finite peroxide treatment and source-backed ammonia/carbonate/cyanide pH.
//!
//! Run with --help for the two explicit models. `treat` closes finite reagent,
//! atom, charge and standard-enthalpy accounts for prescribed extents. `speciate`
//! solves analytical pools at one common ideal pH and a fixed solvent basis;
//! it does not reinterpret that approximation as a finite-solvent H/O ledger.

use std::{error::Error, path::Path};

use fs_matdb::{NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy};
use fs_qty::{
    Amount, Dims, MolarMass, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{
    SpeciesId,
    stoichiometry::{
        ConservedReaction, FormationEnthalpySource,
        minerals::MineralReaction,
        network::ReactionNetwork,
        thermodynamics::{
            FormationGibbsSource,
            acid_base::{
                AcidDissociation, AcidFamily, IdealAcidBaseSystem, SpeciationOptions, SpectatorIon,
                WaterAutoprotolysis,
            },
        },
    },
};

const REFERENCE: &str = "source-elemental-reference-nbs1989";
const USAGE: &str = "usage:
  cyanide_treatment treat PACK_DIR CN_MOL H2O2_MOL WATER_MOL OX_EXTENT_MOL DECOMP_EXTENT_MOL HYDROLYSIS_EXTENT_MOL
  cyanide_treatment speciate PACK_DIR AMMONIA_TOTAL_MOL INORGANIC_CARBON_TOTAL_MOL UNCOMPLEXED_CYANIDE_TOTAL_MOL NREF_MOL SODIUM_MOL CHLORIDE_MOL

treat: all products initially zero; oxidation and decomposition share one stock:
  CN- + H2O2(aq) -> OCN- + H2O(l)
  2 H2O2(aq) -> 2 H2O(l) + O2(g)
A subsequent explicit stage consumes only the available cyanate and water:
  OCN- + 2 H2O(l) -> NH4+ + CO3(2-)
Extents are moles of each complete equation. No rate or treatment efficiency
is inferred; nitrogen remains in cyanide, cyanate or ammonium.

speciate: ammonia means NH4+ + NH3; inorganic carbon means the source's formal
H2CO3 component + HCO3- + CO3(2-); uncomplexed cyanide means HCN + CN-.
All families share one solved pH. Sodium and chloride are explicit spectators.
NREF_MOL is the amount at unit activity, solvent kg times 1 mol/kg for these
NBS sources, not solution volume. Water activity is fixed at one. Cyanate,
peroxide, bound cyanide, gas loss and nonideal activities are outside this model.

Both modes require the phase-specific NBS packs documented in
data/matdb/seed-v1/MINING_CYANIDE_TREATMENT.md at 298.15 K and 100000 Pa.";

struct Source {
    species: SpeciesId,
    pack: NormalizedPack,
    at: QueryPoint,
    molar_mass: f64,
}

fn point(phase: &str, ion: bool, carbonic: bool) -> Result<QueryPoint, Box<dyn Error>> {
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
    if carbonic {
        flags.push(("source-carbonic-component-co2-plus-water", 1.0));
    }
    for (axis, value) in flags {
        at = at.with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)?;
    }
    Ok(at)
}

fn load(directory: &Path, species: &SpeciesId) -> Result<Source, Box<dyn Error>> {
    let (slug, phase, ion) = match species.as_str() {
        "H2O2(aq)" => (
            "hydrogen-peroxide-aqueous-nbs1989",
            "source-phase-aqueous",
            false,
        ),
        "OCN[-](aq)" => ("cyanate-ion-aqueous-nbs1989", "source-phase-aqueous", true),
        "NH3(aq)" => ("ammonia-aqueous-nbs1989", "source-phase-aqueous", false),
        "NH4[+](aq)" => ("ammonium-ion-aqueous-nbs1989", "source-phase-aqueous", true),
        "CN[-](aq)" => ("cyanide-ion-aqueous-nbs1989", "source-phase-aqueous", true),
        "HCN(aq)" => (
            "hydrogen-cyanide-aqueous-nbs1989",
            "source-phase-aqueous",
            false,
        ),
        "H+(aq)" => ("hydrogen-ion-aqueous-nbs1989", "source-phase-aqueous", true),
        "OH[-](aq)" => (
            "hydroxide-ion-aqueous-nbs1989",
            "source-phase-aqueous",
            true,
        ),
        "CO3[2-](aq)" => (
            "carbonate-ion-aqueous-nbs1989",
            "source-phase-aqueous",
            true,
        ),
        "HCO3[-](aq)" => (
            "bicarbonate-ion-aqueous-nbs1989",
            "source-phase-aqueous",
            true,
        ),
        "H2CO3(aq)" => (
            "carbonic-component-aqueous-nbs1989",
            "source-phase-aqueous",
            false,
        ),
        "O2(g)" => (
            "oxygen-gas-nbs1989-thermochemistry",
            "source-phase-gas",
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
    let decoded = NormalizedPack::from_bytes(&bytes)?;
    let pack = NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes)?;
    if pack.pack_id() != slug {
        return Err(format!("{} must contain pack {slug}", path.display()).into());
    }
    let at = point(phase, ion, species.as_str() == "H2CO3(aq)")?;
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

fn heat_sources(sources: &[Source]) -> Vec<FormationEnthalpySource<'_>> {
    sources
        .iter()
        .map(|s| FormationEnthalpySource {
            species: &s.species,
            pack: &s.pack,
            at: &s.at,
            reference_axis: REFERENCE,
        })
        .collect()
}

fn gibbs_sources<'a>(
    sources: &'a [Source],
    reaction: &ConservedReaction,
) -> Result<Vec<FormationGibbsSource<'a>>, Box<dyn Error>> {
    reaction
        .terms()
        .iter()
        .map(|term| {
            let s = sources
                .iter()
                .find(|s| s.species == term.species)
                .ok_or("missing explicit species source")?;
            Ok(FormationGibbsSource {
                species: &s.species,
                pack: &s.pack,
                at: &s.at,
                reference_axis: REFERENCE,
            })
        })
        .collect()
}

fn receipts(sources: &[Source]) {
    println!("source_species,pack_id,content_hash");
    for source in sources {
        println!(
            "{},{},{:?}",
            source.species.as_str(),
            source.pack.pack_id(),
            source.pack.content_hash()
        );
    }
}

fn treat(directory: &Path, values: &[f64]) -> Result<(), Box<dyn Error>> {
    let reactions = [
        MineralReaction::CyanidePeroxideOxidation,
        MineralReaction::HydrogenPeroxideDecomposition,
        MineralReaction::CyanateHydrolysisAmmoniumCarbonate,
    ]
    .into_iter()
    .map(|recipe| recipe.build())
    .collect::<Result<Vec<_>, _>>()?;
    let stage_one: Vec<_> = reactions
        .iter()
        .zip([values[3], values[4], 0.0])
        .map(|(r, n)| (r.id().clone(), Amount::new(n)))
        .collect();
    let stage_two: Vec<_> = reactions
        .iter()
        .zip([0.0, 0.0, values[5]])
        .map(|(r, n)| (r.id().clone(), Amount::new(n)))
        .collect();
    let network = ReactionNetwork::new(reactions)?;
    let stocks = [
        ("CN[-](aq)", values[0]),
        ("H2O2(aq)", values[1]),
        ("H2O(l)", values[2]),
    ];
    let initial: Vec<_> = network
        .species()
        .iter()
        .map(|species| {
            let n = stocks
                .iter()
                .find(|(name, _)| *name == species.as_str())
                .map_or(0.0, |(_, n)| *n);
            (species.clone(), Amount::new(n))
        })
        .collect();
    let sources = network
        .species()
        .iter()
        .map(|s| load(directory, s))
        .collect::<Result<Vec<_>, _>>()?;
    let heat = heat_sources(&sources);
    let oxidized = network.advance_with_standard_enthalpy(&initial, &stage_one, &heat)?;
    let hydrolyzed = network.advance_with_standard_enthalpy(
        &oxidized.inventory.final_inventory,
        &stage_two,
        &heat,
    )?;

    // Resolve both stages and every diagnostic before publishing accepted output.
    let mut diagnostics = Vec::new();
    for reaction in network.reactions() {
        let active: Vec<_> = heat
            .iter()
            .filter(|source| {
                reaction
                    .terms()
                    .iter()
                    .any(|term| &term.species == source.species)
            })
            .map(|source| FormationEnthalpySource {
                species: source.species,
                pack: source.pack,
                at: source.at,
                reference_axis: source.reference_axis,
            })
            .collect();
        let caloric = reaction.standard_equilibrium(&active)?;
        let direct = reaction.standard_gibbs_change(&gibbs_sources(&sources, reaction)?)?;
        diagnostics.push((
            reaction.id().as_str(),
            caloric.enthalpy().molar().value(),
            caloric.molar_entropy().value(),
            caloric.molar_gibbs().value(),
            direct.molar().value(),
        ));
    }
    println!("Prescribed finite treatment stages; NBS 298.15 K, 100000 Pa.");
    println!(
        "reaction,standard_delta_h_J_per_mol,standard_delta_s_J_per_mol_K,standard_delta_g_from_h_s_J_per_mol,direct_standard_delta_g_J_per_mol"
    );
    for (id, h, s, derived_g, direct_g) in diagnostics {
        println!("{id},{h},{s},{derived_g},{direct_g}");
    }
    println!("stage,reaction,extent_mol,standard_delta_h_J");
    for (name, stage) in [
        ("oxidation-and-decomposition", &oxidized),
        ("hydrolysis", &hydrolyzed),
    ] {
        for reaction in &stage.reactions {
            println!(
                "{name},{},{},{}",
                reaction.reaction.as_str(),
                reaction.extent.value(),
                reaction.enthalpy.value()
            );
        }
    }
    println!("species,initial_mol,after_oxidation_mol,after_hydrolysis_mol,final_species_mass_kg");
    let mut initial_mass = 0.0;
    let mut final_mass = 0.0;
    for (index, source) in sources.iter().enumerate() {
        let before = initial[index].1.value();
        let after = hydrolyzed.inventory.final_inventory[index].1.value();
        initial_mass += before * source.molar_mass;
        final_mass += after * source.molar_mass;
        println!(
            "{},{},{},{},{}",
            source.species.as_str(),
            before,
            oxidized.inventory.final_inventory[index].1.value(),
            after,
            after * source.molar_mass
        );
    }
    println!(
        "Total standard chemical enthalpy change [J]: {}",
        oxidized.enthalpy.value() + hydrolyzed.enthalpy.value()
    );
    println!(
        "Source-mass residual, final minus initial [kg]: {}",
        final_mass - initial_mass
    );
    receipts(&sources);
    println!(
        "Nitrogen is retained in the named species. These prescribed extents supply no treatment rate, WAD/total-cyanide removal, residual-water quality, catalyst demand or operating heat duty."
    );
    Ok(())
}

fn speciate(directory: &Path, values: &[f64]) -> Result<(), Box<dyn Error>> {
    let id = |name: &str| SpeciesId::new(name).map_err(Box::<dyn Error>::from);
    let sources = [
        "H2O(l)",
        "H+(aq)",
        "OH[-](aq)",
        "NH4[+](aq)",
        "NH3(aq)",
        "H2CO3(aq)",
        "HCO3[-](aq)",
        "CO3[2-](aq)",
        "HCN(aq)",
        "CN[-](aq)",
    ]
    .into_iter()
    .map(|name| load(directory, &id(name)?))
    .collect::<Result<Vec<_>, Box<dyn Error>>>()?;
    let water_reaction = MineralReaction::WaterDissociation.build()?;
    let water_sources = gibbs_sources(&sources, &water_reaction)?;
    let water = WaterAutoprotolysis::from_sources(
        water_reaction,
        id("H2O(l)")?,
        id("H+(aq)")?,
        id("OH[-](aq)")?,
        &water_sources,
    )?;
    let step = |recipe: MineralReaction,
                acid: &str,
                base: &str|
     -> Result<AcidDissociation, Box<dyn Error>> {
        let reaction = recipe.build()?;
        let evidence = gibbs_sources(&sources, &reaction)?;
        Ok(AcidDissociation::from_sources(
            reaction,
            id(acid)?,
            id(base)?,
            id("H+(aq)")?,
            &evidence,
        )?)
    };
    let system = IdealAcidBaseSystem::new(
        water,
        vec![
            AcidFamily::new(
                "ammonia",
                vec![step(
                    MineralReaction::AmmoniumDissociation,
                    "NH4[+](aq)",
                    "NH3(aq)",
                )?],
            )?,
            AcidFamily::new(
                "inorganic-carbon",
                vec![
                    step(
                        MineralReaction::CarbonicComponentDissociation,
                        "H2CO3(aq)",
                        "HCO3[-](aq)",
                    )?,
                    step(
                        MineralReaction::BicarbonateDissociation,
                        "HCO3[-](aq)",
                        "CO3[2-](aq)",
                    )?,
                ],
            )?,
            AcidFamily::new(
                "uncomplexed-cyanide",
                vec![step(
                    MineralReaction::HydrogenCyanideDissociation,
                    "HCN(aq)",
                    "CN[-](aq)",
                )?],
            )?,
        ],
    )?;
    let spectators = [("Na[+](aq)", 1, values[4]), ("Cl[-](aq)", -1, values[5])]
        .into_iter()
        .map(|(name, charge, n)| {
            Ok(SpectatorIon {
                species: id(name)?,
                charge,
                amount: Amount::new(n),
            })
        })
        .collect::<Result<Vec<_>, Box<dyn Error>>>()?;
    let answer = system.speciate(
        &[
            ("ammonia".into(), Amount::new(values[0])),
            ("inorganic-carbon".into(), Amount::new(values[1])),
            ("uncomplexed-cyanide".into(), Amount::new(values[2])),
        ],
        &spectators,
        Amount::new(values[3]),
        SpeciationOptions::default(),
        || false,
    )?;
    println!(
        "Common ideal pH at NBS 298.15 K, 100000 Pa: {}",
        answer.ph.value()
    );
    println!("family,species,amount_mol,fraction,relative_activity");
    for family in &answer.families {
        for form in &family.forms {
            println!(
                "{},{},{},{},{}",
                family.id,
                form.species.as_str(),
                form.amount.value(),
                form.fraction.value(),
                form.activity.value()
            );
        }
    }
    println!("Proton amount [mol]: {}", answer.proton_amount.value());
    println!(
        "Hydroxide amount [mol]: {}",
        answer.hydroxide_amount.value()
    );
    println!("Spectator species,amount_mol,charge_number");
    for spectator in &answer.spectators {
        println!(
            "{},{},{}",
            spectator.species.as_str(),
            spectator.amount.value(),
            spectator.charge
        );
    }
    println!(
        "Charge residual [mol charge equivalents]: {}",
        answer.charge_balance_residual.value()
    );
    println!(
        "Maximum absolute ln(Q/K): {}",
        answer.max_log_mass_action_residual.value()
    );
    receipts(&sources);
    println!(
        "Fixed solvent basis and unit water activity; analytical totals and explicit counterions are caller inputs. Cyanate/peroxide chemistry, finite-solvent H/O balance, gas loss and nonideal corrections are not represented by this speciation result."
    );
    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() == 1 && matches!(args[0].as_str(), "--help" | "-h") {
        println!("{USAGE}");
        return Ok(());
    }
    if args.len() != 8 {
        return Err(USAGE.into());
    }
    let values = args[2..]
        .iter()
        .map(|text| text.parse::<f64>())
        .collect::<Result<Vec<_>, _>>()?;
    if values
        .iter()
        .any(|value| !value.is_finite() || *value < 0.0)
    {
        return Err("every input amount and extent must be finite and nonnegative".into());
    }
    match args[0].as_str() {
        "treat" => treat(Path::new(&args[1]), &values),
        "speciate" => speciate(Path::new(&args[1]), &values),
        _ => Err(USAGE.into()),
    }
}
