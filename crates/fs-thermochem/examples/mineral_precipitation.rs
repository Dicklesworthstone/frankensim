//! Simultaneous iron/goethite partition and aqueous pH from original USGS data.
//!
//! Compile the five profiles listed by --help into PACK_DIR, then run:
//! ```text
//! cargo run -p fs-thermochem --example mineral_precipitation -- \
//!   --materials PACK_DIR 0.001 0.003 0 10
//! ```
//! The four quantities are total Fe component, fixed chloride, fixed sodium
//! and amount at unit aqueous activity, all in mol. The example's explicit
//! inputs resemble a dilute ferric-chloride analytical inventory; they are
//! authored model scenarios, not a reagent recipe or measured plant solution.
//! Raising chloride at fixed Fe models separately declared acid equivalents;
//! raising sodium models separately declared base equivalents. Salt dissolution,
//! commercial reagent purity, kinetics and dosing are not calculated.
//!
//! Only free Fe3+, water ions, fixed spectator ions and the selected goethite
//! phase participate. Solvent basis and water/solid unit activities stay fixed.
//! Other iron complexes, competing solids and nonideal activities are omitted.
#![deny(unsafe_code)]

use fs_matdb::{NormalizedPack, QueryPoint};
use fs_qty::{
    Amount, Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{
    SpeciesId,
    stoichiometry::{
        minerals::MineralReaction,
        thermodynamics::{
            FormationGibbsSource,
            acid_base::{
                IdealAcidBaseSystem, SpectatorIon, WaterAutoprotolysis,
                precipitation::{HydrolysisPrecipitation, PrecipitationOptions},
            },
        },
    },
};
use std::{error::Error, path::Path};

const REFERENCE: &str = "source-elemental-reference-usgs-b2131";
const PACKS: [&str; 5] = [
    "iron-iii-ion-aqueous-usgs-b2131",
    "goethite-solid-usgs-b2131-thermochemistry",
    "hydrogen-ion-aqueous-usgs-b2131",
    "hydroxide-ion-aqueous-usgs-b2131",
    "water-liquid-usgs-b2131",
];
const SPECIES: [&str; 5] = ["Fe[3+](aq)", "FeO(OH)(s)", "H+(aq)", "OH[-](aq)", "H2O(l)"];
const USAGE: &str = "usage: mineral_precipitation --materials PACK_DIR TOTAL_FE_MOL CHLORIDE_MOL SODIUM_MOL AQUEOUS_UNIT_ACTIVITY_MOL";

fn point(index: usize) -> Result<QueryPoint, Box<dyn Error>> {
    let mut at = QueryPoint::new()
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
    let mut flags = vec![
        (REFERENCE, 1.0),
        ("source-standard-state", 1.0),
        ("engineering-reference-model", 1.0),
        ("source-specimen-qualified", 0.0),
    ];
    match index {
        1 => flags.extend([
            ("source-phase-solid", 1.0),
            ("source-phase-goethite", 1.0),
            ("source-pure-species", 1.0),
        ]),
        4 => flags.extend([("source-phase-liquid", 1.0), ("source-pure-species", 1.0)]),
        _ => flags.extend([
            ("source-phase-aqueous", 1.0),
            ("source-aqueous-unit-activity-infinite-dilution", 1.0),
            ("source-individual-ion", 1.0),
            ("source-ionic-reference-hydrogen-zero", 1.0),
        ]),
    }
    for (axis, value) in flags {
        at = at.with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)?;
    }
    Ok(at)
}

fn load(directory: &Path) -> Result<HydrolysisPrecipitation, Box<dyn Error>> {
    let packs = PACKS
        .iter()
        .map(|name| {
            let bytes = std::fs::read(directory.join(format!("{name}.fsmatpk")))?;
            let decoded = NormalizedPack::from_bytes(&bytes)?;
            let pack = NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes)?;
            if pack.pack_id() != *name {
                return Err(format!(
                    "expected source profile {name}, received {}",
                    pack.pack_id()
                )
                .into());
            }
            Ok(pack)
        })
        .collect::<Result<Vec<_>, Box<dyn Error>>>()?;
    let points = (0..5).map(point).collect::<Result<Vec<_>, _>>()?;
    let ids = SPECIES
        .iter()
        .map(|s| SpeciesId::new(*s))
        .collect::<Result<Vec<_>, _>>()?;
    let sources = |indices: &[usize]| {
        indices
            .iter()
            .map(|&i| FormationGibbsSource {
                species: &ids[i],
                pack: &packs[i],
                at: &points[i],
                reference_axis: REFERENCE,
            })
            .collect::<Vec<_>>()
    };
    let water = WaterAutoprotolysis::from_sources(
        MineralReaction::WaterDissociation.build()?,
        ids[4].clone(),
        ids[2].clone(),
        ids[3].clone(),
        &sources(&[4, 2, 3]),
    )?;
    let aqueous = IdealAcidBaseSystem::new(water, vec![])?;
    Ok(HydrolysisPrecipitation::from_sources(
        aqueous,
        MineralReaction::FerricHydrolysisGoethite.build()?,
        ids[0].clone(),
        ids[1].clone(),
        &sources(&[0, 4, 1, 2]),
    )?)
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    if args.is_empty() || args.as_slice() == ["--help"] {
        println!("{USAGE}");
        println!("Compile these source profiles as PACK_DIR/PACK_ID.fsmatpk:");
        for pack in PACKS {
            println!("  {pack}");
        }
        println!(
            "Use: cargo run -p xtask -- matdb-pack --manifest data/matdb/seed-v1/PACK_ID/manifest.tsv --out PACK_DIR/PACK_ID.fsmatpk"
        );
        println!(
            "All stock arguments are finite nonnegative moles; the unit-activity amount is positive."
        );
        println!("TOTAL_FE_MOL is free Fe3+ plus the selected goethite on a one-Fe formula basis.");
        println!(
            "The fixed solvent basis gives a=n/n_ref, gamma=1; pure water and goethite have unit activity."
        );
        println!(
            "Chloride and sodium are explicit fixed spectator stocks, not inferred counterions or doses."
        );
        return Ok(());
    }
    if args.len() != 6 || args[0] != "--materials" {
        return Err(USAGE.into());
    }
    let model = load(Path::new(&args[1]))?;
    let numbers = args[2..]
        .iter()
        .map(|s| s.parse::<f64>())
        .collect::<Result<Vec<_>, _>>()?;
    let total = Amount::new(numbers[0]);
    let spectators = [("Cl[-](aq)", -1, numbers[1]), ("Na[+](aq)", 1, numbers[2])]
        .into_iter()
        .map(|(name, charge, amount)| {
            Ok(SpectatorIon {
                species: SpeciesId::new(name)?,
                charge,
                amount: Amount::new(amount),
            })
        })
        .collect::<Result<Vec<_>, Box<dyn Error>>>()?;
    let basis = Amount::new(numbers[3]);
    let answer = model.speciate(
        total,
        &[],
        &spectators,
        basis,
        PrecipitationOptions::default(),
        || false,
    )?;
    println!("Coupled Fe3+(aq) + 2 H2O(l) <=> goethite + 3 H+(aq) and water autoprotolysis.");
    println!(
        "Original USGS reference: 298.15 K, 100000 Pa; fixed ideal aqueous basis {} mol.",
        basis.value()
    );
    println!("All aqueous activities use a=n/n_ref and gamma=1; water and goethite activity=1.");
    println!("Total Fe component [mol]: {}", total.value());
    println!(
        "Fixed chloride/sodium [mol]: {} / {}",
        numbers[1], numbers[2]
    );
    println!("Status: {:?}", answer.status);
    println!("pH=-log10(aH): {:.12}", answer.aqueous.ph.value());
    println!("Free Fe3+ [mol]: {:.16e}", answer.dissolved_amount.value());
    println!("Goethite [mol]: {:.16e}", answer.solid_amount.value());
    println!(
        "Free H+ / OH- [mol]: {:.12e} / {:.12e}",
        answer.aqueous.proton_amount.value(),
        answer.aqueous.hydroxide_amount.value()
    );
    println!(
        "Component balance residual [mol]: {:.3e}",
        answer.component_balance_residual.value()
    );
    println!(
        "Aqueous charge residual [mol equivalents]: {:.3e}",
        answer.aqueous.charge_balance_residual.value()
    );
    match answer.log_saturation_ratio {
        Some(value) => println!("Natural log of dissolution IAP/K: {:.12e}", value.value()),
        None => println!("Phase saturation is undefined for an empty Fe component."),
    }
    if let Some(spread) = answer.log_saturation_numerical_spread {
        println!(
            "Remaining log-saturation variation over the inner pH bracket: {:.12e}",
            spread.value()
        );
    }
    println!(
        "Numerical dissolved bracket [mol]: [{:.16e}, {:.16e}]",
        answer.dissolved_bracket[0].value(),
        answer.dissolved_bracket[1].value()
    );
    println!(
        "Outer bisections / aqueous solves / inner pH bisections: {} / {} / {}",
        answer.iterations, answer.aqueous_solves, answer.aqueous_iterations
    );
    println!(
        "An absent phase retains its undersaturation inequality; no mineral mass is invented."
    );
    println!(
        "No other Fe complexes, competing minerals, nonideal activities, kinetics or finite-solvent H/O balance are included."
    );
    for pack in PACKS {
        println!("Original material source: {pack}");
    }
    Ok(())
}
