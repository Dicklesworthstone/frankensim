//! Source-bound finite calcite/gypsum competition with one common pH and Ca stock.
//! Analytical totals include their solid contents; the solvent basis is explicit.
//! See data/matdb/seed-v1/MINING_LIMESTONE_SULFATE.md for sources and model limits.
#![deny(unsafe_code)]

use fs_matdb::{NormalizedPack, QueryPoint};
use fs_qty::{
    Amount, Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{
    SpeciesId,
    stoichiometry::{
        ConservedReaction,
        minerals::MineralReaction,
        thermodynamics::{
            FormationGibbsSource,
            acid_base::{
                AcidDissociation, AcidFamily, IdealAcidBaseSystem, SpectatorIon,
                WaterAutoprotolysis,
                precipitation::salt::{
                    SaltPrecipitation,
                    assemblage::{
                        CommonCationSaltSystem, SaltAssemblageOptions, SaltAssemblageSpeciation,
                    },
                },
            },
        },
    },
};
use std::{error::Error, path::Path};

pub(crate) const REFERENCE: &str = "source-elemental-reference-nbs1989";
pub(crate) const CARBON: &str = "inorganic-carbon";
pub(crate) const SULFATE: &str = "sulfate";
pub(crate) const CALCITE: &str = "CaCO3(s)";
pub(crate) const GYPSUM: &str = "CaSO4.2H2O(s)";
pub(crate) const ROWS: [(&str, &str, &str); 11] = [
    (CALCITE, "calcite-solid-nbs1989-thermochemistry", "calcite"),
    (
        GYPSUM,
        "gypsum-selenite-solid-nbs1989-thermochemistry",
        "gypsum",
    ),
    ("Ca[2+](aq)", "calcium-ion-aqueous-nbs1989", "ion"),
    (
        "H2CO3(aq)",
        "carbonic-component-aqueous-nbs1989",
        "carbonic",
    ),
    ("HCO3[-](aq)", "bicarbonate-ion-aqueous-nbs1989", "ion"),
    ("CO3[2-](aq)", "carbonate-ion-aqueous-nbs1989", "ion"),
    ("HSO4[-](aq)", "hydrogen-sulfate-ion-aqueous-nbs1989", "ion"),
    ("SO4[2-](aq)", "sulfate-ion-aqueous-nbs1989", "ion"),
    ("H2O(l)", "water-liquid-nbs1989-thermochemistry", "liquid"),
    ("H+(aq)", "hydrogen-ion-aqueous-nbs1989", "ion"),
    ("OH[-](aq)", "hydroxide-ion-aqueous-nbs1989", "ion"),
];
const USAGE: &str = "usage: limestone_sulfate_equilibrium PACK_DIR [CA_TOTAL_MOL C_TOTAL_MOL S_TOTAL_MOL NA_MOL CL_MOL UNIT_ACTIVITY_AMOUNT_MOL]";

pub(crate) struct Source {
    pub(crate) id: SpeciesId,
    pub(crate) pack: NormalizedPack,
    pub(crate) at: QueryPoint,
}

fn point(phase: &str) -> Result<QueryPoint, Box<dyn Error>> {
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
        "ion" | "carbonic" => {
            flags.extend([
                ("source-phase-aqueous", 1.0),
                ("source-aqueous-unit-activity-infinite-dilution", 1.0),
            ]);
            if phase == "ion" {
                flags.extend([
                    ("source-individual-ion", 1.0),
                    ("source-ionic-reference-hydrogen-zero", 1.0),
                ]);
            } else {
                flags.push(("source-carbonic-component-co2-plus-water", 1.0));
            }
        }
        "liquid" => flags.extend([("source-phase-liquid", 1.0), ("source-pure-species", 1.0)]),
        "calcite" | "gypsum" => {
            flags.extend([("source-phase-solid", 1.0), ("source-pure-species", 1.0)]);
            if phase == "calcite" {
                flags.push(("source-phase-calcite", 1.0));
            } else {
                flags.extend([
                    ("source-phase-gypsum", 1.0),
                    ("source-variety-selenite", 1.0),
                ]);
            }
        }
        _ => return Err("unknown original source phase".into()),
    }
    for (name, value) in flags {
        at = at.with_quantity(name, QuantitySpec::dimensional(Dims::NONE), value)?;
    }
    Ok(at)
}

pub(crate) fn load_sources(directory: &Path) -> Result<Vec<Source>, Box<dyn Error>> {
    ROWS.iter()
        .map(|(id, name, phase)| {
            let bytes = std::fs::read(directory.join(format!("{name}.fsmatpk")))?;
            let raw = NormalizedPack::from_bytes(&bytes)?;
            let pack = NormalizedPack::from_bytes_verified(raw.content_hash(), &bytes)?;
            if pack.pack_id() != *name {
                return Err(format!("expected {name}, found {}", pack.pack_id()).into());
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

pub(crate) fn model_from_sources(
    sources: &[Source],
) -> Result<CommonCationSaltSystem, Box<dyn Error>> {
    let id = SpeciesId::new;
    let r = MineralReaction::WaterDissociation.build()?;
    let water = WaterAutoprotolysis::from_sources(
        r.clone(),
        id("H2O(l)")?,
        id("H+(aq)")?,
        id("OH[-](aq)")?,
        &selected(&r, sources),
    )?;
    let mut carbon = Vec::new();
    for (recipe, acid, base) in [
        (
            MineralReaction::CarbonicComponentDissociation,
            "H2CO3(aq)",
            "HCO3[-](aq)",
        ),
        (
            MineralReaction::BicarbonateDissociation,
            "HCO3[-](aq)",
            "CO3[2-](aq)",
        ),
    ] {
        let r = recipe.build()?;
        carbon.push(AcidDissociation::from_sources(
            r.clone(),
            id(acid)?,
            id(base)?,
            id("H+(aq)")?,
            &selected(&r, sources),
        )?);
    }
    let r = MineralReaction::HydrogenSulfateDissociation.build()?;
    let sulfate = AcidDissociation::from_sources(
        r.clone(),
        id("HSO4[-](aq)")?,
        id("SO4[2-](aq)")?,
        id("H+(aq)")?,
        &selected(&r, sources),
    )?;
    let aqueous = IdealAcidBaseSystem::new(
        water,
        vec![
            AcidFamily::new(CARBON, carbon)?,
            AcidFamily::new(SULFATE, vec![sulfate])?,
        ],
    )?;
    let mut phases = Vec::new();
    for (recipe, solid, family) in [
        (MineralReaction::CalcitePrecipitation, CALCITE, CARBON),
        (MineralReaction::GypsumPrecipitation, GYPSUM, SULFATE),
    ] {
        let r = recipe.build()?;
        phases.push(SaltPrecipitation::from_sources(
            aqueous.clone(),
            r.clone(),
            id("Ca[2+](aq)")?,
            id(solid)?,
            family,
            &selected(&r, sources),
        )?);
    }
    Ok(CommonCationSaltSystem::new(phases)?)
}

/// Analytical component stocks; calcium/carbon/sulfur INCLUDE any selected solids.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Charge {
    pub(crate) calcium: Amount,
    pub(crate) carbon: Amount,
    pub(crate) sulfur: Amount,
    pub(crate) sodium: Amount,
    pub(crate) chloride: Amount,
    pub(crate) basis: Amount,
}
impl Charge {
    pub(crate) fn families(self) -> [(String, Amount); 2] {
        [(CARBON.into(), self.carbon), (SULFATE.into(), self.sulfur)]
    }
    pub(crate) fn spectators(self) -> Result<[SpectatorIon; 2], Box<dyn Error>> {
        Ok([
            SpectatorIon {
                species: SpeciesId::new("Na[+](aq)")?,
                charge: 1,
                amount: self.sodium,
            },
            SpectatorIon {
                species: SpeciesId::new("Cl[-](aq)")?,
                charge: -1,
                amount: self.chloride,
            },
        ])
    }
}
pub(crate) fn demonstration_charges() -> [Charge; 2] {
    let initial = Charge {
        calcium: Amount::new(0.05),
        carbon: Amount::new(0.03),
        sulfur: Amount::new(0.05),
        sodium: Amount::new(0.10),
        chloride: Amount::new(0.06),
        basis: Amount::new(1.0),
    };
    [
        initial,
        Charge {
            chloride: Amount::new(0.065),
            ..initial
        },
    ]
}
pub(crate) fn solve(
    model: &CommonCationSaltSystem,
    charge: Charge,
) -> Result<SaltAssemblageSpeciation, Box<dyn Error>> {
    Ok(model.speciate(
        charge.calcium,
        &charge.families(),
        &charge.spectators()?,
        charge.basis,
        SaltAssemblageOptions::default(),
        || false,
    )?)
}

fn display(label: &str, a: &SaltAssemblageSpeciation) {
    println!(
        "case={label} pH={:.10} dissolved_Ca_mol={:.12e}",
        a.aqueous.ph.value(),
        a.dissolved_cation.value()
    );
    for family in &a.aqueous.families {
        for form in &family.forms {
            println!(
                "aqueous={} amount_mol={:.12e} activity={:.12e}",
                form.species.as_str(),
                form.amount.value(),
                form.activity.value()
            );
        }
    }
    for phase in &a.phases {
        println!(
            "solid={} amount_mol={:.12e} status={:?} ln_IAP_over_K={:?}",
            phase.solid.as_str(),
            phase.amount.value(),
            phase.status,
            phase.log_saturation_ratio.map(|v| v.value())
        );
    }
    println!(
        "calcium_balance_mol={:.3e} charge_balance_mol={:.3e}",
        a.cation_balance_residual.value(),
        a.aqueous.charge_balance_residual.value()
    );
    for family in &a.families {
        println!(
            "family={} dissolved_mol={:.12e} in_solid_mol={:.12e} balance_mol={:.3e}",
            family.family,
            family.dissolved.value(),
            family.in_solid.value(),
            family.residual.value()
        );
    }
    println!(
        "sweeps={} coordinate_solves={} aqueous_solves={} pH_iterations={}",
        a.sweeps, a.coordinate_solves, a.aqueous_solves, a.aqueous_iterations
    );
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() != 1 && args.len() != 7 {
        return Err(USAGE.into());
    }
    let sources = load_sources(Path::new(&args[0]))?;
    let model = model_from_sources(&sources)?;
    let cases = if args.len() == 1 {
        demonstration_charges()
            .into_iter()
            .zip(["initial", "after_0.005_mol_HCl"])
            .collect::<Vec<_>>()
    } else {
        let values = args[1..]
            .iter()
            .map(|s| s.parse::<f64>())
            .collect::<Result<Vec<_>, _>>()?;
        vec![(
            Charge {
                calcium: Amount::new(values[0]),
                carbon: Amount::new(values[1]),
                sulfur: Amount::new(values[2]),
                sodium: Amount::new(values[3]),
                chloride: Amount::new(values[4]),
                basis: Amount::new(values[5]),
            },
            "custom",
        )]
    };
    // Admit and solve EVERY requested case before emitting a successful report.
    let results = cases
        .into_iter()
        .map(|(charge, label)| Ok((label, solve(&model, charge)?)))
        .collect::<Result<Vec<_>, Box<dyn Error>>>()?;
    println!(
        "model=ideal_fixed_solvent_calcite_gypsum temperature_K=298.15 pressure_Pa=100000 source_uncertainty=Unstated"
    );
    for (label, answer) in results {
        display(label, &answer);
    }
    Ok(())
}
