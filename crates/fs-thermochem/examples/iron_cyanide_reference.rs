//! Explicitly selected legacy NBS iron-cyanide reference calculation.
//!
//! All four iron profiles are original source data, not fitted plant chemistry.
//! The NBS Fe(II) association constant materially disagrees with later critical
//! literature. This executable therefore requires --reference nbs1989-legacy.
//! Two independently declared Fe oxidation-state pools compete for one finite
//! HCN/CN-/bound-CN inventory. No redox, hydrolysis, precipitation, photolysis,
//! protonated complexes, solid-ore access or kinetic rates are inferred.
//! See data/matdb/seed-v1/MINING_IRON_CYANIDE_REFERENCE.md for source limitations.
#![deny(unsafe_code)]

use std::{error::Error, path::Path};

use fs_matdb::{NormalizedPack, QueryPoint};
use fs_qty::{
    Amount, Dimensionless, Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{
    ChargeVector, ElementId, ElementalMatrix, ReactionId, SpeciesId, StoichiometricMatrix,
    stoichiometry::{
        ConservedReaction,
        minerals::MineralReaction,
        thermodynamics::{
            FormationGibbsSource,
            acid_base::{
                AcidDissociation, AcidFamily, IdealAcidBaseSystem, SpectatorIon,
                WaterAutoprotolysis,
                complexation::{
                    IdealLigandSystem, LigandAssociation, LigandSpeciation,
                    LigandSpeciationOptions, MetalLigandFamily,
                },
            },
        },
    },
};

pub(crate) const REFERENCE: &str = "source-elemental-reference-nbs1989";
pub(crate) const LIGAND: &str = "uncomplexed-cyanide";
pub(crate) const FE_II: &str = "declared-iron-ii";
pub(crate) const FE_III: &str = "declared-iron-iii";
pub(crate) const PACKS: [&str; 9] = [
    "hydrogen-cyanide-aqueous-nbs1989",
    "hydrogen-ion-aqueous-nbs1989",
    "cyanide-ion-aqueous-nbs1989",
    "hydroxide-ion-aqueous-nbs1989",
    "water-liquid-nbs1989-thermochemistry",
    "iron-ii-ion-aqueous-nbs1989",
    "iron-iii-ion-aqueous-nbs1989",
    "hexacyanoferrate-ii-ion-aqueous-nbs1989",
    "hexacyanoferrate-iii-ion-aqueous-nbs1989",
];
pub(crate) const SPECIES: [&str; 9] = [
    "HCN(aq)",
    "H+(aq)",
    "CN[-](aq)",
    "OH[-](aq)",
    "H2O(l)",
    "Fe[2+](aq)",
    "Fe[3+](aq)",
    "Fe(CN)6[4-](aq)",
    "Fe(CN)6[3-](aq)",
];
const USAGE: &str = "usage: iron_cyanide_reference --reference nbs1989-legacy --materials PACK_DIR --amounts TOTAL_CN_MOL FE_II_MOL FE_III_MOL NA_MOL CL_MOL UNIT_ACTIVITY_AMOUNT_MOL";

pub(crate) struct Source {
    pub(crate) id: SpeciesId,
    pub(crate) pack: NormalizedPack,
    pub(crate) at: QueryPoint,
}

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
    if index == 4 {
        flags.extend([("source-phase-liquid", 1.0), ("source-pure-species", 1.0)]);
    } else {
        flags.extend([
            ("source-phase-aqueous", 1.0),
            ("source-aqueous-unit-activity-infinite-dilution", 1.0),
        ]);
        if index != 0 {
            flags.extend([
                ("source-individual-ion", 1.0),
                ("source-ionic-reference-hydrogen-zero", 1.0),
            ]);
        }
    }
    for (name, value) in flags {
        at = at.with_quantity(name, QuantitySpec::dimensional(Dims::NONE), value)?;
    }
    Ok(at)
}

pub(crate) fn load_sources(directory: &Path) -> Result<Vec<Source>, Box<dyn Error>> {
    PACKS
        .iter()
        .enumerate()
        .map(|(index, name)| {
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
            Ok(Source {
                id: SpeciesId::new(SPECIES[index])?,
                pack,
                at: point(index)?,
            })
        })
        .collect()
}

pub(crate) fn selected<'a>(
    reaction: &ConservedReaction,
    sources: &'a [Source],
) -> Vec<FormationGibbsSource<'a>> {
    reaction
        .terms()
        .iter()
        .filter_map(|term| {
            sources
                .iter()
                .find(|source| source.id == term.species)
                .map(|source| FormationGibbsSource {
                    species: &source.id,
                    pack: &source.pack,
                    at: &source.at,
                    reference_axis: REFERENCE,
                })
        })
        .collect()
}

pub(crate) fn association(oxidation: u8) -> Result<ConservedReaction, Box<dyn Error>> {
    let (parent, complex) = match oxidation {
        2 => (5, 7),
        3 => (6, 8),
        _ => return Err("iron oxidation-state pool must be 2 or 3".into()),
    };
    let species = [SPECIES[parent], SPECIES[2], SPECIES[complex]]
        .into_iter()
        .map(SpeciesId::new)
        .collect::<Result<Vec<_>, _>>()?;
    let elements = ["Fe", "C", "N"]
        .into_iter()
        .map(ElementId::new)
        .collect::<Result<Vec<_>, _>>()?;
    let composition = ElementalMatrix::new(
        elements,
        species.clone(),
        vec![vec![1, 0, 1], vec![0, 1, 6], vec![0, 1, 6]],
    )?;
    let reaction = ReactionId::new(format!("iron-{oxidation}-hexacyanide-association"))?;
    let coefficients = StoichiometricMatrix::new(
        species.clone(),
        vec![reaction.clone()],
        vec![vec![-1], vec![-6], vec![1]],
    )?;
    let charges = ChargeVector::new(
        species,
        vec![i128::from(oxidation), -1, i128::from(oxidation) - 6],
    )?;
    Ok(ConservedReaction::new(
        &composition,
        &coefficients,
        &charges,
        reaction,
    )?)
}

pub(crate) fn model_from_sources(sources: &[Source]) -> Result<IdealLigandSystem, Box<dyn Error>> {
    let water_reaction = MineralReaction::WaterDissociation.build()?;
    let water = WaterAutoprotolysis::from_sources(
        water_reaction.clone(),
        SpeciesId::new(SPECIES[4])?,
        SpeciesId::new(SPECIES[1])?,
        SpeciesId::new(SPECIES[3])?,
        &selected(&water_reaction, sources),
    )?;
    let acid_reaction = MineralReaction::HydrogenCyanideDissociation.build()?;
    let acid = AcidDissociation::from_sources(
        acid_reaction.clone(),
        SpeciesId::new(SPECIES[0])?,
        SpeciesId::new(SPECIES[2])?,
        SpeciesId::new(SPECIES[1])?,
        &selected(&acid_reaction, sources),
    )?;
    let aqueous = IdealAcidBaseSystem::new(water, vec![AcidFamily::new(LIGAND, vec![acid])?])?;
    let mut families = Vec::new();
    for (oxidation, family, parent, complex) in [(2, FE_II, 5, 7), (3, FE_III, 6, 8)] {
        let reaction = association(oxidation)?;
        let step = LigandAssociation::from_sources(
            reaction.clone(),
            SpeciesId::new(SPECIES[parent])?,
            SpeciesId::new(SPECIES[2])?,
            SpeciesId::new(SPECIES[complex])?,
            &selected(&reaction, sources),
        )?;
        families.push(MetalLigandFamily::new(family, vec![step])?);
    }
    Ok(IdealLigandSystem::new(aqueous, LIGAND, families)?)
}

/// Explicit input amounts: cyanide, Fe(II), Fe(III), Na+, Cl-, amount at unit activity.
/// Na/Cl are independent caller stocks; no counterions are manufactured.
pub(crate) fn solve(
    model: &IdealLigandSystem,
    amounts: [f64; 6],
) -> Result<LigandSpeciation, Box<dyn Error>> {
    if amounts.iter().any(|x| !x.is_finite() || *x < 0.0) || amounts[5] <= 0.0 {
        return Err("all amounts must be finite and nonnegative; activity-reference amount must be positive".into());
    }
    let [cyanide, fe_ii, fe_iii, sodium, chloride, basis] = amounts;
    Ok(model.speciate(
        Amount::new(cyanide),
        &[
            (FE_II.into(), Amount::new(fe_ii)),
            (FE_III.into(), Amount::new(fe_iii)),
        ],
        &[],
        &[
            SpectatorIon {
                species: SpeciesId::new("Na[+](aq)")?,
                charge: 1,
                amount: Amount::new(sodium),
            },
            SpectatorIon {
                species: SpeciesId::new("Cl[-](aq)")?,
                charge: -1,
                amount: Amount::new(chloride),
            },
        ],
        Amount::new(basis),
        LigandSpeciationOptions {
            log_ligand_tolerance: Dimensionless::new(1e-12),
            absolute_ligand_tolerance: Amount::new(1e-17),
            relative_ligand_tolerance: Dimensionless::new(1e-13),
            ..LigandSpeciationOptions::default()
        },
        || false,
    )?)
}

pub(crate) fn parse_amounts(args: &[String]) -> Result<(&Path, [f64; 6]), Box<dyn Error>> {
    if args.len() != 11
        || args[0] != "--reference"
        || args[1] != "nbs1989-legacy"
        || args[2] != "--materials"
        || args[4] != "--amounts"
    {
        return Err(USAGE.into());
    }
    let mut amounts = [0.0; 6];
    for (value, argument) in amounts.iter_mut().zip(&args[5..]) {
        *value = argument.parse::<f64>()?;
    }
    if amounts.iter().any(|x| !x.is_finite() || *x < 0.0) || amounts[5] <= 0.0 {
        return Err("all amounts must be finite and nonnegative; activity-reference amount must be positive".into());
    }
    Ok((Path::new(&args[3]), amounts))
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    if args.as_slice() == ["--help"] {
        println!("{USAGE}");
        println!("All amounts are mol. Activities are a=n/n_ref, gamma=1; water activity is one.");
        println!("Compile each profile as PACK_DIR/PACK_ID.fsmatpk:");
        for pack in PACKS {
            println!("  {pack}");
        }
        println!(
            "NBS1989 is an explicitly selected legacy reference, not validated capture chemistry."
        );
        return Ok(());
    }
    let (directory, amounts) = parse_amounts(&args)?;
    let sources = load_sources(directory)?;
    let model = model_from_sources(&sources)?;
    let solution = solve(&model, amounts)?;
    println!("REFERENCE=nbs1989-legacy; T=298.15 K; p=100000 Pa; ideal aqueous activities.");
    println!(
        "CAUTION: NBS Fe(II) log10(beta)=45.623610; Keefe/Miller (1996) use thermodynamic log10(beta)=35.397940."
    );
    println!(
        "Source disagreement is unresolved. This is a reference calculation, not a validated capture prediction."
    );
    println!(
        "Fe(II)/Fe(III) pools are fixed declarations; no redox, hydrolysis, precipitation, photolysis or kinetics."
    );
    println!(
        "n_ref={:.12e} mol; each Na+/Cl- spectator stock is supplied explicitly.",
        amounts[5]
    );
    for family in model.families() {
        println!(
            "{} NBS log10(beta)={:.12}",
            family.id(),
            family.steps()[0]
                .standard()
                .log_equilibrium_constant()
                .value()
                / std::f64::consts::LN_10
        );
    }
    println!(
        "NBS_REFERENCE pH={:.12}; free_CN_mol={:.12e}; HCN_plus_free_CN_mol={:.12e}; bound_CN_mol={:.12e}",
        solution.aqueous.ph.value(),
        solution.free_ligand.value(),
        solution.uncomplexed_ligand.value(),
        solution.bound_ligand.value()
    );
    for family in &solution.metals {
        for member in &family.members {
            println!(
                "NBS_REFERENCE {} amount_mol={:.12e}; charge={}; cyanides_per_iron={}",
                member.species.as_str(),
                member.amount.value(),
                member.charge,
                member.ligand_count
            );
        }
        println!(
            "NBS_REFERENCE {} iron_balance_residual_mol={:.12e}",
            family.id,
            family.members.iter().map(|m| m.amount.value()).sum::<f64>() - family.total.value()
        );
    }
    println!(
        "NBS_REFERENCE cyanide_balance_residual_mol={:.12e}; charge_balance_residual_mol={:.12e}; max_log_mass_action_residual={:.12e}",
        solution.ligand_balance_residual.value(),
        solution.aqueous.charge_balance_residual.value(),
        solution.max_log_mass_action_residual.value()
    );
    Ok(())
}
