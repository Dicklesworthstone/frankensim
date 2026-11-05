//! Source-bound Cu2S phase partition, explicit separation and cyanide reuse.
//! The three equilibrium stages share finite Cu/S/CN stocks; no kinetic,
//! gas-loss, filtration-efficiency or plant-performance calibration is implied.
//! See data/matdb/seed-v1/MINING_CYANIDE_RECOVERY.md.
#![deny(unsafe_code)]

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
                AcidDissociation, AcidFamily, IdealAcidBaseSystem, SpeciationOptions, SpectatorIon,
                WaterAutoprotolysis,
                complexation::{
                    IdealLigandSystem, LigandAssociation, LigandSpeciationOptions,
                    MetalLigandFamily,
                    precipitation::{
                        CopperSulfideOptions, CopperSulfidePrecipitation, CopperSulfideSpeciation,
                    },
                },
            },
        },
    },
};
use std::{error::Error, path::Path};

pub(crate) const REFERENCE: &str = "source-elemental-reference-nbs1989";
pub(crate) const COPPER: &str = "dissolved-copper-i";
pub(crate) const SULFIDE: &str = "dissolved-hydrogen-sulfide";
pub(crate) const CYANIDE: &str = "uncomplexed-cyanide";
pub(crate) const SOLID: &str = "Cu2S(s-nbs-alpha)";
pub(crate) const ROWS: [(&str, &str, &str); 12] = [
    ("H2S(aq)", "hydrogen-sulfide-aqueous-nbs1989", "aqueous"),
    ("HS[-](aq)", "hydrosulfide-ion-aqueous-nbs1989", "ion"),
    (SOLID, "copper-i-sulfide-alpha-solid-nbs1989", "solid"),
    ("HCN(aq)", "hydrogen-cyanide-aqueous-nbs1989", "aqueous"),
    ("H+(aq)", "hydrogen-ion-aqueous-nbs1989", "ion"),
    ("CN[-](aq)", "cyanide-ion-aqueous-nbs1989", "ion"),
    ("OH[-](aq)", "hydroxide-ion-aqueous-nbs1989", "ion"),
    ("H2O(l)", "water-liquid-nbs1989-thermochemistry", "liquid"),
    ("Cu[+](aq)", "copper-i-ion-aqueous-nbs1989", "ion"),
    (
        "Cu(CN)2[-](aq)",
        "dicyanocuprate-i-ion-aqueous-nbs1989",
        "ion",
    ),
    (
        "Cu(CN)3[2-](aq)",
        "tricyanocuprate-i-ion-aqueous-nbs1989",
        "ion",
    ),
    (
        "Cu(CN)4[3-](aq)",
        "tetracyanocuprate-i-ion-aqueous-nbs1989",
        "ion",
    ),
];
const USAGE: &str = "usage: copper_cyanide_recovery PACK_DIR CN_MOL CU_MOL FEED_NA_MOL FEED_CL_MOL NAHS_MOL HCL_MOL NAOH_MOL UNIT_ACTIVITY_AMOUNT_MOL";

pub(crate) struct Source {
    pub(crate) id: SpeciesId,
    pub(crate) pack: NormalizedPack,
    pub(crate) at: QueryPoint,
}
pub(crate) fn point(phase: &str) -> Result<QueryPoint, Box<dyn Error>> {
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
    match phase {
        "ion" | "aqueous" => {
            flags.extend([
                ("source-phase-aqueous", 1.0),
                ("source-aqueous-unit-activity-infinite-dilution", 1.0),
            ]);
            if phase == "ion" {
                flags.extend([
                    ("source-individual-ion", 1.0),
                    ("source-ionic-reference-hydrogen-zero", 1.0),
                ]);
            }
        }
        "liquid" => flags.extend([("source-phase-liquid", 1.0), ("source-pure-species", 1.0)]),
        "solid" => flags.extend([
            ("source-phase-solid", 1.0),
            ("source-pure-species", 1.0),
            ("source-phase-copper-i-sulfide-nbs-alpha", 1.0),
        ]),
        _ => return Err("unknown source phase".into()),
    }
    for (name, value) in flags {
        at = at.with_quantity(name, QuantitySpec::dimensional(Dims::NONE), value)?;
    }
    Ok(at)
}
pub(crate) fn load_sources(directory: &Path) -> Result<Vec<Source>, Box<dyn Error>> {
    ROWS.iter()
        .map(|(id, pack_id, phase)| {
            let bytes = std::fs::read(directory.join(format!("{pack_id}.fsmatpk")))?;
            let raw = NormalizedPack::from_bytes(&bytes)?;
            let pack = NormalizedPack::from_bytes_verified(raw.content_hash(), &bytes)?;
            if pack.pack_id() != *pack_id {
                return Err(format!("expected {pack_id}, found {}", pack.pack_id()).into());
            }
            Ok(Source {
                id: SpeciesId::new(*id)?,
                pack,
                at: point(phase)?,
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
                .find(|s| s.id == term.species)
                .map(|s| FormationGibbsSource {
                    species: &s.id,
                    pack: &s.pack,
                    at: &s.at,
                    reference_axis: REFERENCE,
                })
        })
        .collect()
}
pub(crate) fn copper_association(level: u8) -> Result<ConservedReaction, Box<dyn Error>> {
    let (parent, product, prior, added) = match level {
        2 => ("Cu[+](aq)", "Cu(CN)2[-](aq)", 0, 2),
        3 => ("Cu(CN)2[-](aq)", "Cu(CN)3[2-](aq)", 2, 1),
        4 => ("Cu(CN)3[2-](aq)", "Cu(CN)4[3-](aq)", 3, 1),
        _ => return Err("copper association level must be 2,3 or4".into()),
    };
    let species = [parent, "CN[-](aq)", product]
        .into_iter()
        .map(SpeciesId::new)
        .collect::<Result<Vec<_>, _>>()?;
    let elements = ["Cu", "C", "N"]
        .into_iter()
        .map(ElementId::new)
        .collect::<Result<Vec<_>, _>>()?;
    let composition = ElementalMatrix::new(
        elements,
        species.clone(),
        vec![
            vec![1, 0, 1],
            vec![prior, 1, prior + added],
            vec![prior, 1, prior + added],
        ],
    )?;
    let reaction = ReactionId::new(format!("copper-cyanide-association-{level}"))?;
    let stoichiometry = StoichiometricMatrix::new(
        species.clone(),
        vec![reaction.clone()],
        vec![vec![-1], vec![-added], vec![1]],
    )?;
    let charges = ChargeVector::new(species, vec![1 - prior, -1, 1 - prior - added])?;
    Ok(ConservedReaction::new(
        &composition,
        &stoichiometry,
        &charges,
        reaction,
    )?)
}
pub(crate) fn model_from_sources(
    sources: &[Source],
) -> Result<CopperSulfidePrecipitation, Box<dyn Error>> {
    let id = SpeciesId::new;
    let water = MineralReaction::WaterDissociation.build()?;
    let water = WaterAutoprotolysis::from_sources(
        water.clone(),
        id("H2O(l)")?,
        id("H+(aq)")?,
        id("OH[-](aq)")?,
        &selected(&water, sources),
    )?;
    let mut families = Vec::new();
    for (name, recipe, acid, base) in [
        (
            CYANIDE,
            MineralReaction::HydrogenCyanideDissociation,
            "HCN(aq)",
            "CN[-](aq)",
        ),
        (
            SULFIDE,
            MineralReaction::HydrogenSulfideDissociation,
            "H2S(aq)",
            "HS[-](aq)",
        ),
    ] {
        let r = recipe.build()?;
        let step = AcidDissociation::from_sources(
            r.clone(),
            id(acid)?,
            id(base)?,
            id("H+(aq)")?,
            &selected(&r, sources),
        )?;
        families.push(AcidFamily::new(name, vec![step])?);
    }
    let aqueous = IdealAcidBaseSystem::new(water, families)?;
    let mut steps = Vec::new();
    for (level, parent, product) in [
        (2, "Cu[+](aq)", "Cu(CN)2[-](aq)"),
        (3, "Cu(CN)2[-](aq)", "Cu(CN)3[2-](aq)"),
        (4, "Cu(CN)3[2-](aq)", "Cu(CN)4[3-](aq)"),
    ] {
        let r = copper_association(level)?;
        steps.push(LigandAssociation::from_sources(
            r.clone(),
            id(parent)?,
            id("CN[-](aq)")?,
            id(product)?,
            &selected(&r, sources),
        )?);
    }
    let ligand = IdealLigandSystem::new(
        aqueous,
        CYANIDE,
        vec![MetalLigandFamily::new(COPPER, steps)?],
    )?;
    let phase = MineralReaction::CopperSulfideHydrosulfidePrecipitation.build()?;
    Ok(CopperSulfidePrecipitation::from_sources(
        ligand,
        phase.clone(),
        COPPER,
        SULFIDE,
        id(SOLID)?,
        &selected(&phase, sources),
    )?)
}

/// All six stocks are independent molar amounts: CN, Cu, sulfide, Na, Cl, n_ref.
pub(crate) fn solve(
    model: &CopperSulfidePrecipitation,
    amounts: [f64; 6],
) -> Result<CopperSulfideSpeciation, Box<dyn Error>> {
    let [cn, cu, s, na, cl, basis] = amounts;
    if amounts.iter().any(|x| !x.is_finite() || *x < 0.0) || basis <= 0.0 {
        return Err(
            "finite nonnegative stocks and positive activity-reference amount required".into(),
        );
    }
    Ok(model.speciate(
        Amount::new(cn),
        &[(COPPER.into(), Amount::new(cu))],
        &[(SULFIDE.into(), Amount::new(s))],
        &[
            SpectatorIon {
                species: SpeciesId::new("Na[+](aq)")?,
                charge: 1,
                amount: Amount::new(na),
            },
            SpectatorIon {
                species: SpeciesId::new("Cl[-](aq)")?,
                charge: -1,
                amount: Amount::new(cl),
            },
        ],
        Amount::new(basis),
        CopperSulfideOptions {
            // Resolve the complete-sequence CN balance more tightly than its
            // reporting tolerance, including nearly complete Cu removal.
            ligand: LigandSpeciationOptions {
                aqueous: SpeciationOptions {
                    ph_tolerance: Dimensionless::new(1e-13),
                    absolute_charge_tolerance: Amount::new(1e-17),
                    relative_charge_tolerance: Dimensionless::new(1e-13),
                    ..SpeciationOptions::default()
                },
                log_ligand_tolerance: Dimensionless::new(1e-11),
                absolute_ligand_tolerance: Amount::new(1e-16),
                relative_ligand_tolerance: Dimensionless::new(1e-11),
                ..LigandSpeciationOptions::default()
            },
            ..CopperSulfideOptions::default()
        },
        || false,
    )?)
}
pub(crate) struct Stages {
    pub(crate) feed: CopperSulfideSpeciation,
    pub(crate) acidified: CopperSulfideSpeciation,
    pub(crate) recycled: CopperSulfideSpeciation,
    pub(crate) separated_solid: f64,
}
/// Amounts: total CN, Cu, feed Na, feed Cl, NaHS addition, HCl addition,
/// NaOH addition, n_ref. Additions are declared ion equivalents at fixed solvent
/// basis, not a dissolution/heat/volume model for commercial reagent products.
pub(crate) fn run_stages(
    model: &CopperSulfidePrecipitation,
    amounts: [f64; 8],
) -> Result<Stages, Box<dyn Error>> {
    if amounts.iter().any(|x| !x.is_finite() || *x < 0.0) || amounts[7] <= 0.0 {
        return Err(
            "finite nonnegative stocks and positive activity-reference amount required".into(),
        );
    }
    let [cn, cu, na, cl, nahs, hcl, naoh, basis] = amounts;
    let feed = solve(model, [cn, cu, 0.0, na, cl, basis])?;
    let acidified = solve(model, [cn, cu, nahs, na + nahs, cl + hcl, basis])?;
    // Explicit ideal removal of the calculated solid is a caller-declared
    // separation action, not a fitted or predicted filtration efficiency.
    let separated_solid = acidified.solid_amount.value();
    let recycled = solve(
        model,
        [
            cn,
            acidified.dissolved_copper.value(),
            acidified.dissolved_sulfide.value(),
            na + nahs + naoh,
            cl + hcl,
            basis,
        ],
    )?;
    // Independent whole-sequence element and ligand checks include both solids.
    let all_solid = separated_solid + recycled.solid_amount.value();
    for (actual, initial) in [
        (recycled.dissolved_copper.value() + 2.0 * all_solid, cu),
        (recycled.dissolved_sulfide.value() + all_solid, nahs),
        (
            recycled.ligand.uncomplexed_ligand.value() + recycled.ligand.bound_ligand.value(),
            cn,
        ),
    ] {
        if !actual.is_finite() || (actual - initial).abs() > 1e-10 * initial.abs().max(1e-12) {
            return Err("whole-sequence Cu/S/CN inventory does not close".into());
        }
    }
    Ok(Stages {
        feed,
        acidified,
        recycled,
        separated_solid,
    })
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() != 9 {
        return Err(USAGE.into());
    }
    let sources = load_sources(Path::new(&args[0]))?;
    let model = model_from_sources(&sources)?;
    let parsed = args[1..]
        .iter()
        .map(|s| s.parse::<f64>())
        .collect::<Result<Vec<_>, _>>()?;
    let amounts: [f64; 8] = parsed.try_into().map_err(|_| USAGE)?;
    let results = run_stages(&model, amounts)?;
    // Publish only after all three stages and whole-sequence balances succeed.
    println!(
        "Selected NBS Cu2S(alpha)/Cu(I)-cyanide/H2S-HS ideal assemblage,298.15K,1bar; admitted pH2..10."
    );
    println!(
        "Declared additions: NaHS={}mol,HCl={}mol,NaOH={}mol; n_ref={}mol (not liquid volume).",
        amounts[4], amounts[5], amounts[6], amounts[7]
    );
    println!(
        "Phase delta_G={}J/mol,ln_K_precip={}",
        model.standard().molar().value(),
        model.standard().log_equilibrium_constant().value()
    );
    for (name, state) in [
        ("feed", &results.feed),
        ("acidified", &results.acidified),
        ("after solid separation and NaOH", &results.recycled),
    ] {
        println!(
            "{name}: pH={:.12},Cu_diss={:.12e}mol,S_diss={:.12e}mol,Cu2S={:.12e}mol,free_CN={:.12e}mol,HCN_plus_CN={:.12e}mol,bound_CN={:.12e}mol,status={:?}",
            state.ligand.aqueous.ph.value(),
            state.dissolved_copper.value(),
            state.dissolved_sulfide.value(),
            state.solid_amount.value(),
            state.ligand.free_ligand.value(),
            state.ligand.uncomplexed_ligand.value(),
            state.ligand.bound_ligand.value(),
            state.status
        );
        println!(
            "  residuals: Cu={:.4e}mol,S={:.4e}mol,CN={:.4e}mol,charge={:.4e}mol; ln_saturation={:?},numerical_spread={:?}",
            state.copper_balance_residual.value(),
            state.sulfide_balance_residual.value(),
            state.ligand.ligand_balance_residual.value(),
            state.ligand.aqueous.charge_balance_residual.value(),
            state.log_saturation_ratio.map(|x| x.value()),
            state.log_saturation_numerical_spread.map(|x| x.value())
        );
    }
    println!(
        "Explicitly separated Cu2S={}mol; any secondary solid remains in the final inventory.",
        results.separated_solid
    );
    println!(
        "No S(2-), other solids/Cu-sulfide complexes, Cu redox, gas loss, kinetic rates or plant calibration. Fixed solvent; no complete H/O ledger."
    );
    for source in sources {
        println!(
            "source {} pack={} hash={:?}",
            source.id.as_str(),
            source.pack.pack_id(),
            source.pack.content_hash()
        );
    }
    Ok(())
}
