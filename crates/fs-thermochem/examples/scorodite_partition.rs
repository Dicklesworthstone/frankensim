//! Finite Fe/As(V) partition with one selected scorodite phase and shared pH.
//!
//! Uses the consistent NMK2014/CODATA source profiles listed by --help. Both
//! analytical component totals include any initially present scorodite. Sodium
//! and chloride are explicit fixed spectators; no counterions are inferred.
//! Activities use a common fixed solvent basis, not solution volume. This
//! selected-phase ideal equilibrium omits ferric hydrolysis/arsenate complexes,
//! competing solids, nonideal activities, adsorption and precipitation kinetics.
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
                AcidDissociation, AcidFamily, IdealAcidBaseSystem, SpectatorIon,
                WaterAutoprotolysis,
                precipitation::{PrecipitationOptions, salt::SaltPrecipitation},
            },
        },
    },
};

const REFERENCE: &str = "source-elemental-reference-nmk2014";
const FAMILY: &str = "uncomplexed-arsenate";
const PACKS: [&str; 9] = [
    "arsenic-acid-aqueous-nmk2014",
    "dihydrogen-arsenate-ion-aqueous-nmk2014",
    "hydrogen-arsenate-ion-aqueous-nmk2014",
    "arsenate-ion-aqueous-nmk2014",
    "iron-iii-ion-aqueous-nmk2014",
    "scorodite-solid-nmk2014",
    "water-liquid-nmk2014",
    "hydrogen-ion-aqueous-nmk2014",
    "hydroxide-ion-aqueous-nmk2014",
];
const SPECIES: [&str; 9] = [
    "H3AsO4(aq)",
    "H2AsO4[-](aq)",
    "HAsO4[2-](aq)",
    "AsO4[3-](aq)",
    "Fe[3+](aq)",
    "FeAsO4.2H2O(s)",
    "H2O(l)",
    "H+(aq)",
    "OH[-](aq)",
];
const USAGE: &str = "usage: scorodite_partition --materials PACK_DIR [--stocks FE_MOL AS_MOL NA_MOL CL_MOL BASIS_MOL]";

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
        5 => flags.extend([
            ("source-phase-solid", 1.0),
            ("source-pure-species", 1.0),
            ("source-phase-scorodite", 1.0),
        ]),
        6 => flags.extend([("source-phase-liquid", 1.0), ("source-pure-species", 1.0)]),
        _ => {
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
    }
    if index < 4 || index == 5 {
        flags.push(("source-arsenic-oxidation-state-v", 1.0));
    }
    for (axis, value) in flags {
        at = at.with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)?;
    }
    Ok(at)
}

fn model(directory: &Path) -> Result<SaltPrecipitation, Box<dyn Error>> {
    let packs = PACKS
        .iter()
        .map(|name| {
            let bytes = std::fs::read(directory.join(format!("{name}.fsmatpk")))?;
            let decoded = NormalizedPack::from_bytes(&bytes)?;
            let pack = NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes)?;
            if pack.pack_id() != *name {
                return Err(format!("expected {name}, received {}", pack.pack_id()).into());
            }
            Ok(pack)
        })
        .collect::<Result<Vec<_>, Box<dyn Error>>>()?;
    let points = (0..PACKS.len()).map(point).collect::<Result<Vec<_>, _>>()?;
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
        ids[6].clone(),
        ids[7].clone(),
        ids[8].clone(),
        &sources(&[6, 7, 8]),
    )?;
    let steps = [
        MineralReaction::ArsenicAcidDissociation,
        MineralReaction::DihydrogenArsenateDissociation,
        MineralReaction::HydrogenArsenateDissociation,
    ]
    .into_iter()
    .enumerate()
    .map(|(i, recipe)| {
        AcidDissociation::from_sources(
            recipe.build()?,
            ids[i].clone(),
            ids[i + 1].clone(),
            ids[7].clone(),
            &sources(&[i, 7, i + 1]),
        )
    })
    .collect::<Result<Vec<_>, _>>()?;
    let aqueous = IdealAcidBaseSystem::new(water, vec![AcidFamily::new(FAMILY, steps)?])?;
    Ok(SaltPrecipitation::from_sources(
        aqueous,
        MineralReaction::ScoroditePrecipitation.build()?,
        ids[4].clone(),
        ids[5].clone(),
        FAMILY,
        &sources(&[4, 3, 6, 5]),
    )?)
}

fn run(model: &SaltPrecipitation, name: &str, stocks: [f64; 5]) -> Result<(), Box<dyn Error>> {
    let [iron, arsenic, sodium, chloride, basis] = stocks;
    let result = model.speciate(
        Amount::new(iron),
        &[(FAMILY.to_owned(), Amount::new(arsenic))],
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
        PrecipitationOptions::default(),
        || false,
    )?;
    println!("scenario: {name}");
    println!("  declared Fe/As/Na/Cl/basis [mol]: {stocks:?}");
    println!(
        "  phase status: {:?}; pH: {:.12}",
        result.status,
        result.aqueous.ph.value()
    );
    println!("  scorodite [mol]: {:.12e}", result.solid_amount.value());
    println!(
        "  dissolved Fe3+ / analytical As(V) [mol]: {:.12e} / {:.12e}",
        result.dissolved_cation.value(),
        result.dissolved_family.value()
    );
    for form in &result.aqueous.families[0].forms {
        println!(
            "  {}: {:.12e} mol; fraction {:.12e}",
            form.species.as_str(),
            form.amount.value(),
            form.fraction.value()
        );
    }
    println!(
        "  Fe/As component residuals [mol]: {:.6e} / {:.6e}",
        result.cation_balance_residual.value(),
        result.family_balance_residual.value()
    );
    println!(
        "  charge residual [mol charge equivalents]: {:.6e}",
        result.aqueous.charge_balance_residual.value()
    );
    println!(
        "  ln(IAP/Ksp): {:?}; numerical pH spread: {:?}",
        result.log_saturation_ratio.map(|v| v.value()),
        result.log_saturation_numerical_spread.map(|v| v.value())
    );
    println!(
        "  limiting dissolved bracket [mol]: {:?}; outer iterations: {}; aqueous solves: {}",
        result.dissolved_limiting_bracket.map(|v| v.value()),
        result.iterations,
        result.aqueous_solves
    );
    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    if args.len() == 1 && args[0] == "--help" {
        println!(
            "{USAGE}\nAll five stocks are explicit nonnegative amounts; BASIS_MOL must be positive."
        );
        println!(
            "Without --stocks, run six declared illustrative scenarios on a 10 kg solvent basis."
        );
        println!("Compile these source profiles as PACK_ID.fsmatpk in PACK_DIR:");
        for pack in PACKS {
            println!("  {pack}");
        }
        return Ok(());
    }
    if !matches!(args.len(), 2 | 8)
        || args[0] != "--materials"
        || (args.len() == 8 && args[2] != "--stocks")
    {
        return Err(USAGE.into());
    }
    let model = model(Path::new(&args[1]))?;
    println!(
        "Selected-phase ideal Fe/As(V)/scorodite equilibrium; NMK2014/CODATA, 298.15 K, 100000 Pa."
    );
    println!(
        "Activities use n/n_ref with fixed unit water activity. No complexation, competing solids, nonideality, rates or plant removal claim."
    );
    if args.len() == 8 {
        let mut stocks = [0.0; 5];
        for (out, arg) in stocks.iter_mut().zip(&args[3..]) {
            *out = arg.parse()?;
        }
        run(&model, "caller-supplied stocks", stocks)?;
    } else {
        for (name, stocks) in [
            (
                "matched iron and arsenic, no sodium",
                [0.02, 0.02, 0.0, 0.06, 10.0],
            ),
            (
                "one sodium equivalent per arsenic",
                [0.02, 0.02, 0.02, 0.06, 10.0],
            ),
            (
                "two sodium equivalents per arsenic",
                [0.02, 0.02, 0.04, 0.06, 10.0],
            ),
            ("iron-limited feed", [0.005, 0.02, 0.0, 0.015, 10.0]),
            ("arsenic-limited feed", [0.04, 0.02, 0.04, 0.12, 10.0]),
            ("trace undersaturated feed", [1e-12, 1e-12, 0.0, 0.1, 10.0]),
        ] {
            run(&model, name, stocks)?;
        }
    }
    Ok(())
}
