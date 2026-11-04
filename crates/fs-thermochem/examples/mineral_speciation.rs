//! Coupled ideal acid/base and water speciation from matched NBS source packs.
//!
//! Run `mineral_speciation --materials PACK_DIR [--family cyanide|sulfate|arsenate]`;
//! cyanide is the default. Each named profile must be
//! compiled there as PACK_ID.fsmatpk. `--help` lists the required profiles.
//! The example uses a fixed 10 kg solvent basis, n_ref=10 mol at unit molal
//! activity, and explicit Na+/Cl- spectator inventories. Neither molarity,
//! solution volume, commercial reagent purity nor measured plant pH is inferred.
//! Water has fixed unit activity: this is not a finite-solvent H/O balance.
//! The analytical pool is HCN+CN- only, excluding gold/copper complexes and WAD-CN.
//! The sulfate mode instead uses HSO4- + SO4(2-), excluding molecular H2SO4
//! and metal-sulfate complexes. It does not model concentrated sulfuric acid.
//! The arsenate mode retains all four uncomplexed As(V) protonation states.
//! As(III), mineral dissolution, redox, adsorption and fixation are separate.
#![deny(unsafe_code)]

use std::{error::Error, path::Path};

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
                AcidDissociation, AcidFamily, IdealAcidBaseSystem, SpeciationOptions, SpectatorIon,
                WaterAutoprotolysis,
            },
        },
    },
};

const REFERENCE: &str = "source-elemental-reference-nbs1989";
const FAMILY: &str = "uncomplexed-cyanide";
const PACKS: [&str; 5] = [
    "hydrogen-cyanide-aqueous-nbs1989",
    "hydrogen-ion-aqueous-nbs1989",
    "cyanide-ion-aqueous-nbs1989",
    "hydroxide-ion-aqueous-nbs1989",
    "water-liquid-nbs1989-thermochemistry",
];
const SPECIES: [&str; 5] = ["HCN(aq)", "H+(aq)", "CN[-](aq)", "OH[-](aq)", "H2O(l)"];
const SULFATE_PACKS: [&str; 5] = [
    "hydrogen-sulfate-ion-aqueous-nbs1989",
    "hydrogen-ion-aqueous-nbs1989",
    "sulfate-ion-aqueous-nbs1989",
    "hydroxide-ion-aqueous-nbs1989",
    "water-liquid-nbs1989-thermochemistry",
];
const SULFATE_SPECIES: [&str; 5] = [
    "HSO4[-](aq)",
    "H+(aq)",
    "SO4[2-](aq)",
    "OH[-](aq)",
    "H2O(l)",
];
const SULFATE_FAMILY: &str = "uncomplexed-sulfate";
const ARSENATE_PACKS: [&str; 7] = [
    "arsenic-acid-aqueous-nbs1989",
    "hydrogen-ion-aqueous-nbs1989",
    "dihydrogen-arsenate-ion-aqueous-nbs1989",
    "hydroxide-ion-aqueous-nbs1989",
    "water-liquid-nbs1989-thermochemistry",
    "hydrogen-arsenate-ion-aqueous-nbs1989",
    "arsenate-ion-aqueous-nbs1989",
];
const ARSENATE_SPECIES: [&str; 7] = [
    "H3AsO4(aq)",
    "H+(aq)",
    "H2AsO4[-](aq)",
    "OH[-](aq)",
    "H2O(l)",
    "HAsO4[2-](aq)",
    "AsO4[3-](aq)",
];
const USAGE: &str =
    "usage: mineral_speciation --materials PACK_DIR [--family cyanide|sulfate|arsenate]";

#[derive(Clone, Copy, PartialEq, Eq)]
enum Family {
    Cyanide,
    Sulfate,
    Arsenate,
}
impl Family {
    fn packs(self) -> &'static [&'static str] {
        match self {
            Self::Cyanide => &PACKS,
            Self::Sulfate => &SULFATE_PACKS,
            Self::Arsenate => &ARSENATE_PACKS,
        }
    }
    fn species(self) -> &'static [&'static str] {
        match self {
            Self::Cyanide => &SPECIES,
            Self::Sulfate => &SULFATE_SPECIES,
            Self::Arsenate => &ARSENATE_SPECIES,
        }
    }
    fn id(self) -> &'static str {
        match self {
            Self::Cyanide => FAMILY,
            Self::Sulfate => SULFATE_FAMILY,
            Self::Arsenate => "uncomplexed-arsenate",
        }
    }
}

fn point(index: usize, family: Family) -> Result<QueryPoint, Box<dyn Error>> {
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
    if index == 4 {
        flags.extend([("source-phase-liquid", 1.0), ("source-pure-species", 1.0)]);
    } else {
        flags.extend([
            ("source-phase-aqueous", 1.0),
            ("source-aqueous-unit-activity-infinite-dilution", 1.0),
        ]);
        if index != 0 || family == Family::Sulfate {
            flags.extend([
                ("source-individual-ion", 1.0),
                ("source-ionic-reference-hydrogen-zero", 1.0),
            ]);
        }
    }
    if family == Family::Arsenate && [0, 2, 5, 6].contains(&index) {
        flags.push(("source-arsenic-oxidation-state-v", 1.0));
    }
    for (axis, value) in flags {
        at = at.with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)?;
    }
    Ok(at)
}

fn formation_sources<'a>(
    indices: &[usize],
    ids: &'a [SpeciesId],
    packs: &'a [NormalizedPack],
    points: &'a [QueryPoint],
) -> Vec<FormationGibbsSource<'a>> {
    indices
        .iter()
        .map(|&i| FormationGibbsSource {
            species: &ids[i],
            pack: &packs[i],
            at: &points[i],
            reference_axis: REFERENCE,
        })
        .collect()
}

fn load_system(directory: &Path, family: Family) -> Result<IdealAcidBaseSystem, Box<dyn Error>> {
    let pack_names = family.packs();
    let species_names = family.species();
    let packs = pack_names
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
    let points = (0..pack_names.len())
        .map(|i| point(i, family))
        .collect::<Result<Vec<_>, _>>()?;
    let ids = species_names
        .iter()
        .map(|s| SpeciesId::new(*s))
        .collect::<Result<Vec<_>, _>>()?;
    let water = WaterAutoprotolysis::from_sources(
        MineralReaction::WaterDissociation.build()?,
        ids[4].clone(),
        ids[1].clone(),
        ids[3].clone(),
        &formation_sources(&[4, 1, 3], &ids, &packs, &points),
    )?;
    let mut steps = vec![(
        match family {
            Family::Cyanide => MineralReaction::HydrogenCyanideDissociation,
            Family::Sulfate => MineralReaction::HydrogenSulfateDissociation,
            Family::Arsenate => MineralReaction::ArsenicAcidDissociation,
        },
        0,
        2,
    )];
    if family == Family::Arsenate {
        steps.extend([
            (MineralReaction::DihydrogenArsenateDissociation, 2, 5),
            (MineralReaction::HydrogenArsenateDissociation, 5, 6),
        ]);
    }
    let dissociations = steps
        .into_iter()
        .map(|(reaction, acid, base)| {
            AcidDissociation::from_sources(
                reaction.build()?,
                ids[acid].clone(),
                ids[base].clone(),
                ids[1].clone(),
                &formation_sources(&[acid, 1, base], &ids, &packs, &points),
            )
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok(IdealAcidBaseSystem::new(
        water,
        vec![AcidFamily::new(family.id(), dissociations)?],
    )?)
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    if args.is_empty() || args.as_slice() == ["--help"] {
        println!("{USAGE}");
        println!("Cyanide source profiles, compiled as PACK_DIR/PACK_ID.fsmatpk:");
        for pack in PACKS {
            println!("  {pack}");
        }
        println!("Sulfate source profiles for --family sulfate:");
        for pack in SULFATE_PACKS {
            println!("  {pack}");
        }
        println!("Arsenate source profiles for --family arsenate:");
        for pack in ARSENATE_PACKS {
            println!("  {pack}");
        }
        println!(
            "Compile each with: cargo run -p xtask -- matdb-pack --manifest data/matdb/seed-v1/PACK_ID/manifest.tsv --out PACK_DIR/PACK_ID.fsmatpk"
        );
        println!(
            "Runs pure water, the chosen acid/base pool, spectator salts, added base and added acid."
        );
        println!("The fixed 10 kg solvent basis gives n_ref=10 mol; a=n/n_ref and gamma=1.");
        println!("Na+ and Cl- are explicit spectators; each pool excludes metal complexes.");
        println!(
            "The sulfate pool is HSO4- + SO4(2-); molecular H2SO4 and concentrated acid are excluded."
        );
        println!("The arsenate pool contains all four As(V) forms, excluding As(III) and solids.");
        println!(
            "Water has unit activity; no finite water balance, redox, metal complexes or kinetics."
        );
        return Ok(());
    }
    if !(args.len() == 2 || args.len() == 4) || args[0] != "--materials" {
        return Err(USAGE.into());
    }
    let family_mode = if args.len() == 4 {
        if args[2] != "--family" {
            return Err(USAGE.into());
        }
        match args[3].as_str() {
            "cyanide" => Family::Cyanide,
            "sulfate" => Family::Sulfate,
            "arsenate" => Family::Arsenate,
            _ => return Err(USAGE.into()),
        }
    } else {
        Family::Cyanide
    };
    let system = load_system(Path::new(&args[1]), family_mode)?;
    let names = family_mode.species();
    let family_id = family_mode.id();
    let options = SpeciationOptions::default();
    let basis = Amount::new(10.0);
    println!("Ideal coupled aqueous speciation at 298.15 K and 100000 Pa.");
    println!("Source evaluation: NIST NBS tables with 1989 errata; direct formation-G route.");
    println!(
        "Fixed solvent mass: 10 kg; amount at unit molal activity: {} mol.",
        basis.value()
    );
    println!("All aqueous activities use a=n/n_ref, gamma=1; water activity is one.");
    println!(
        "Source-derived pKa({})={:.12}; pKw={:.12}.",
        names[0],
        -system.families()[0].steps()[0]
            .standard()
            .log_equilibrium_constant()
            .value()
            / core::f64::consts::LN_10,
        -system.water().standard().log_equilibrium_constant().value() / core::f64::consts::LN_10,
    );
    println!(
        "Spectator amounts are fixed caller declarations, not salt dissolution or dosing predictions."
    );
    for step in system.families()[0].steps().iter().skip(1) {
        println!(
            "Source-derived pKa({})={:.12}.",
            step.acid().as_str(),
            -step.standard().log_equilibrium_constant().value() / core::f64::consts::LN_10,
        );
    }
    let scenarios = if family_mode == Family::Arsenate {
        vec![
            ("pure water; empty uncomplexed As(V) pool", 0.0, 0.0, 0.0),
            ("arsenic acid; no fixed counterions", 0.1, 0.0, 0.0),
            ("one positive counterion per As(V)", 0.1, 0.1, 0.0),
            ("two positive counterions per As(V)", 0.1, 0.2, 0.0),
            ("three positive counterions per As(V)", 0.1, 0.3, 0.0),
            ("two-counterion pool plus acid equivalents", 0.1, 0.2, 0.05),
        ]
    } else if family_mode == Family::Sulfate {
        vec![
            ("pure water; empty uncomplexed sulfate pool", 0.0, 0.0, 0.0),
            (
                "uncomplexed sulfate acid pool; no fixed counterions",
                0.1,
                0.0,
                0.0,
            ),
            ("bisulfate salt with one Na+ per sulfur", 0.1, 0.1, 0.0),
            ("sulfate salt with two Na+ per sulfur", 0.1, 0.2, 0.0),
            (
                "sulfate salt plus 0.02 mol strong-base equivalents",
                0.1,
                0.22,
                0.0,
            ),
            (
                "sulfate salt plus 0.05 mol strong-acid equivalents",
                0.1,
                0.2,
                0.05,
            ),
        ]
    } else {
        vec![
            ("pure water; empty uncomplexed cyanide pool", 0.0, 0.0, 0.0),
            ("neutral molecular HCN feed", 0.1, 0.0, 0.0),
            ("cyanide salt with matching Na+", 0.1, 0.1, 0.0),
            (
                "cyanide salt plus 0.02 mol strong-base equivalents",
                0.1,
                0.12,
                0.0,
            ),
            (
                "cyanide salt plus 0.05 mol strong-acid equivalents",
                0.1,
                0.1,
                0.05,
            ),
        ]
    };
    for (label, total, sodium, chloride) in scenarios {
        let mut spectators = Vec::new();
        for (name, charge, amount) in [("Na[+](aq)", 1, sodium), ("Cl[-](aq)", -1, chloride)] {
            if amount > 0.0 {
                spectators.push(SpectatorIon {
                    species: SpeciesId::new(name)?,
                    charge,
                    amount: Amount::new(amount),
                });
            }
        }
        let answer = system.speciate(
            &[(family_id.to_owned(), Amount::new(total))],
            &spectators,
            basis,
            options,
            || false,
        )?;
        let family = &answer.families[0];
        let acid = family
            .forms
            .iter()
            .find(|f| f.species.as_str() == names[0])
            .unwrap();
        let base = family
            .forms
            .iter()
            .find(|f| f.species.as_str() == names[2])
            .unwrap();
        let closure = family.forms.iter().map(|f| f.amount.value()).sum::<f64>() - total;
        println!();
        println!("Scenario: {label}");
        println!("Declared pool: {total} mol; Na+: {sodium} mol; Cl-: {chloride} mol.");
        println!("pH=-log10(aH): {:.10}", answer.ph.value());
        println!(
            "{}: {:.12e} mol; {}: {:.12e} mol.",
            names[0],
            acid.amount.value(),
            names[2],
            base.amount.value()
        );
        if total > 0.0 {
            println!(
                "{} fraction of the uncomplexed pool: {:.10}",
                names[2],
                base.fraction.value()
            );
        }
        for form in family
            .forms
            .iter()
            .filter(|form| form.species.as_str() != names[0] && form.species.as_str() != names[2])
        {
            println!(
                "{}: {:.12e} mol; fraction of the uncomplexed pool: {:.10}.",
                form.species.as_str(),
                form.amount.value(),
                form.fraction.value()
            );
        }
        println!(
            "Free H+: {:.12e} mol; free OH-: {:.12e} mol.",
            answer.proton_amount.value(),
            answer.hydroxide_amount.value()
        );
        println!("Sum of family species minus declared pool: {closure:.3e} mol.");
        println!(
            "Charge-balance residual: {:.3e} mol charge equivalents.",
            answer.charge_balance_residual.value()
        );
        println!(
            "Maximum absolute ln(Q/K) residual: {:.3e}.",
            answer.max_log_mass_action_residual.value()
        );
        println!(
            "Numerical pH bracket: [{:.12}, {:.12}]; bisections: {}.",
            answer.ph_bracket[0].value(),
            answer.ph_bracket[1].value(),
            answer.iterations
        );
    }
    println!();
    println!("Each uncomplexed pool excludes metal complexes; these are ideal-model predictions.");
    println!(
        "No gas transfer, volatilization, precipitation, redox or leach-rate model is solved."
    );
    println!(
        "The solvent basis and water activity stay fixed; no closed solvent H/O balance is claimed."
    );
    Ok(())
}
