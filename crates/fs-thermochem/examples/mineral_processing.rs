//! Evaluate lime slaking from the three compiled, source-backed material packs.
//!
//! Usage: `cargo run -p fs-thermochem --example mineral_processing -- \
//! cao.fsmatpk water.fsmatpk caoh2.fsmatpk 56.077`
//! Or explicit malachite inventory accounting, with all four inputs in moles:
//! `cargo run -p fs-thermochem --example mineral_processing -- \
//! --malachite-inventory 10 12 100 3`
//!
//! The final argument is kilograms of pure CaO actually reacted. It is not
//! commercial feed mass; feed purity and conversion require separate evidence.
//! This is standard-state bookkeeping at 298.15 K and 1 bar, not a plant
//! temperature, reaction-rate, water-dosing or recovery prediction.

use std::{error::Error, path::Path};

use fs_matdb::{MatDbError, NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Amount, Dims, Mass, MolarMass, QuantitySpec};
use fs_thermochem::{
    SpeciesId,
    stoichiometry::{FormationEnthalpySource, minerals::MineralReaction},
};

const REFERENCE: &str = "source-elemental-reference-janaf";
const PACK_IDS: [&str; 3] = [
    "calcium-oxide-solid-nist-mining",
    "water-liquid-janaf-mining",
    "calcium-hydroxide-solid-nist-mining",
];

fn load(path: &Path, expected_id: &str) -> Result<NormalizedPack, Box<dyn Error>> {
    let bytes = std::fs::read(path)?;
    let decoded = NormalizedPack::from_bytes(&bytes)?;
    let pack = NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes)?;
    // This example explicitly associates known source profiles with phases;
    // arbitrary property packs do not authenticate a chemical formula.
    if pack.pack_id() != expected_id {
        return Err(format!("expected pack {expected_id}, received {}", pack.pack_id()).into());
    }
    Ok(pack)
}

fn reference_point(liquid: bool) -> Result<QueryPoint, MatDbError> {
    let mut point = QueryPoint::new()
        .with_quantity(
            "temperature",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            298.15,
        )?
        .with_quantity(
            "pressure",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            100_000.0,
        )?;
    for (axis, value) in [
        ("source-standard-state", 1.0),
        (REFERENCE, 1.0),
        ("source-pure-species", 1.0),
        ("engineering-reference-model", 1.0),
        ("source-specimen-qualified", 0.0),
        (
            if liquid {
                "source-phase-liquid"
            } else {
                "source-phase-solid"
            },
            1.0,
        ),
    ] {
        point = point.with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)?;
    }
    Ok(point)
}

fn malachite_inventory(args: &[String]) -> Result<(), Box<dyn Error>> {
    if args.len() != 5 {
        return Err("usage: mineral_processing --malachite-inventory MALACHITE_MOL PROTON_MOL WATER_MOL EXTENT_MOL".into());
    }
    let reaction = MineralReaction::MalachiteAcidLeach.build()?;
    let initial = [
        (
            SpeciesId::new("Cu2CO3(OH)2(s)")?,
            Amount::new(args[1].parse()?),
        ),
        (SpeciesId::new("H+(aq)")?, Amount::new(args[2].parse()?)),
        (SpeciesId::new("H2O(l)")?, Amount::new(args[3].parse()?)),
        (SpeciesId::new("Cu[2+](aq)")?, Amount::new(0.0)),
        (SpeciesId::new("CO2(aq)")?, Amount::new(0.0)),
    ];
    let extent = Amount::new(args[4].parse()?);
    let bound = reaction.limiting_extent(&initial[..2])?;
    let final_inventory = reaction.advance_inventory(&initial, extent)?;
    println!("Cu2CO3(OH)2(s) + 4 H+(aq) -> 2 Cu[2+](aq) + CO2(aq) + 3 H2O(l)");
    println!(
        "material-availability extent bound [mol]: {}",
        bound.value()
    );
    println!("caller-declared extent [mol]: {}", extent.value());
    for (species, remaining) in final_inventory {
        let before = initial
            .iter()
            .find(|(name, _)| *name == species)
            .ok_or("final species is absent from initial inventory")?
            .1;
        println!(
            "{} [mol]: {} -> {}",
            species.as_str(),
            before.value(),
            remaining.value()
        );
    }
    println!(
        "Extent is an input, not predicted recovery. Proton inventory is not a commercial acid dose."
    );
    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args
        .first()
        .is_some_and(|arg| arg == "--malachite-inventory")
    {
        return malachite_inventory(&args);
    }
    if args.len() != 4 {
        return Err(
            "usage: mineral_processing CAO_PACK WATER_PACK CAOH2_PACK REACTED_PURE_CAO_KG".into(),
        );
    }
    let reacted_kg: f64 = args[3].parse()?;
    let packs = [
        load(Path::new(&args[0]), PACK_IDS[0])?,
        load(Path::new(&args[1]), PACK_IDS[1])?,
        load(Path::new(&args[2]), PACK_IDS[2])?,
    ];
    let species = [
        SpeciesId::new("CaO(s)")?,
        SpeciesId::new("H2O(l)")?,
        SpeciesId::new("Ca(OH)2(s)")?,
    ];
    let points = [
        reference_point(false)?,
        reference_point(true)?,
        reference_point(false)?,
    ];
    let sources = std::array::from_fn::<_, 3, _>(|i| FormationEnthalpySource {
        species: &species[i],
        pack: &packs[i],
        at: &points[i],
        reference_axis: REFERENCE,
    });
    let reaction = MineralReaction::LimeHydration.build()?;
    let heat = reaction.standard_enthalpy_change(&sources)?;
    let molar_mass = packs[0].claims().query_typed(
        &PropertyKey::with_quantity(
            "molar-mass",
            QuantitySpec::dimensional(Dims([0, 1, 0, 0, 0, -1])),
        ),
        &points[0],
        SelectionPolicy::SingleClaimOnly,
    )?;
    packs[0].claims().verify_receipt(&molar_mass.receipt)?;
    let extent = reaction.extent_from_reacted_mass(
        &species[0],
        Mass::new(reacted_kg),
        MolarMass::new(molar_mass.evidence.value.value),
    )?;
    let amounts = reaction.amounts(extent)?;
    let water = amounts
        .iter()
        .find(|amount| amount.term.species == species[1])
        .ok_or("water is absent from the slaking equation")?;

    println!("CaO(s) + H2O(l) -> Ca(OH)2(s), 298.15 K, 100000 Pa");
    println!("reacted pure CaO mass [kg]: {reacted_kg}");
    println!(
        "source CaO molar mass [kg/mol]: {}",
        molar_mass.evidence.value.value
    );
    println!("reacted CaO amount [mol]: {}", extent.value());
    println!(
        "stoichiometric water consumed [mol]: {}",
        water.amount.value()
    );
    println!(
        "standard reaction enthalpy [J/mol extent]: {}",
        heat.molar().value()
    );
    println!(
        "standard enthalpy change [J]: {}",
        heat.at_extent(extent)?.value()
    );
    for pack in &packs {
        println!("source pack: {} ({})", pack.pack_id(), pack.content_hash());
    }
    println!("Negative enthalpy is exothermic. Source uncertainty is unstated.");
    println!(
        "This uses actually reacted pure CaO; it does not predict conversion, plant temperature or water dosing."
    );
    Ok(())
}
