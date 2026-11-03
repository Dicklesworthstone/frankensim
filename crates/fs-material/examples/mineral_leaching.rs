//! Run: cargo run -p fs-material --example mineral_leaching
//! Add `--thermochemistry PACK_DIR` to use the five actual USGS dolomite,
//! H+, Ca2+, Mg2+ and bicarbonate packs for reaction heat at the computed extent.
//!
//! Synthetic finite-batch dolomite/acid calculation. The original conserved
//! equation supplies species consumption/production; an explicit authored
//! bilinear availability law supplies finite-time progress. This k is a fixture,
//! not a published dolomite rate, calibrated ore model, or surface mechanism.

use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, NormalizedPack, Provenance,
    QueryPoint,
};
use fs_material::reaction_batch::{
    BILINEAR_BATCH_LAW, BilinearBatchReaction, BilinearRateConstant, ELAPSED_TIME_AXIS,
    RATE_PARAMETER, REACTANT_A_CONCENTRATION_AXIS, REACTANT_B_CONCENTRATION_AXIS, TEMPERATURE_AXIS,
    VOLUME_AXIS,
};
use fs_qty::{
    Amount, AmountConcentration, Dims, QuantitySpec, Temperature, Time, Volume,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{
    SpeciesId,
    stoichiometry::{
        ConservedReaction, FormationEnthalpySource, ReactionEnthalpy, minerals::MineralReaction,
    },
};
use std::{collections::BTreeMap, error::Error, path::Path};

const USAGE: &str = "usage: mineral_leaching [--thermochemistry PACK_DIR]";

fn reference_enthalpy(
    directory: &Path,
    reaction: &ConservedReaction,
) -> Result<ReactionEnthalpy, Box<dyn Error>> {
    const PACKS: [&str; 5] = [
        "dolomite-solid-usgs-b2131-thermochemistry",
        "hydrogen-ion-aqueous-usgs-b2131",
        "calcium-ion-aqueous-usgs-b2131",
        "magnesium-ion-aqueous-usgs-b2131",
        "bicarbonate-ion-aqueous-usgs-b2131",
    ];
    const REFERENCE: &str = "source-elemental-reference-usgs-b2131";
    let mut packs = Vec::new();
    let mut points = Vec::new();
    let ids = [
        "CaMg(CO3)2(s)",
        "H+(aq)",
        "Ca[2+](aq)",
        "Mg[2+](aq)",
        "HCO3[-](aq)",
    ]
    .map(SpeciesId::new)
    .into_iter()
    .collect::<Result<Vec<_>, _>>()?;
    for (i, name) in PACKS.iter().enumerate() {
        let bytes = std::fs::read(directory.join(format!("{name}.fsmatpk")))?;
        let decoded = NormalizedPack::from_bytes(&bytes)?;
        let pack = NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes)?;
        if pack.pack_id() != *name || pack.to_bytes() != bytes {
            return Err(format!("expected canonical source pack {name}").into());
        }
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
                QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::Pressure,
                    ValueForm::Static,
                )),
                100_000.0,
            )?;
        for (axis, value) in [
            ("source-standard-state", 1.0),
            (REFERENCE, 1.0),
            ("engineering-reference-model", 1.0),
            ("source-specimen-qualified", 0.0),
        ] {
            point = point.with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)?;
        }
        let phase: &[&str] = if i == 0 {
            &[
                "source-phase-solid",
                "source-phase-dolomite",
                "source-pure-species",
            ]
        } else {
            &[
                "source-phase-aqueous",
                "source-individual-ion",
                "source-aqueous-unit-activity-infinite-dilution",
                "source-ionic-reference-hydrogen-zero",
            ]
        };
        for axis in phase {
            point = point.with_quantity(*axis, QuantitySpec::dimensional(Dims::NONE), 1.0)?;
        }
        println!(
            "Thermochemical source: {} ({})",
            pack.pack_id(),
            pack.content_hash()
        );
        packs.push(pack);
        points.push(point);
    }
    let sources = (0..PACKS.len())
        .map(|i| FormationEnthalpySource {
            species: &ids[i],
            pack: &packs[i],
            at: &points[i],
            reference_axis: REFERENCE,
        })
        .collect::<Vec<_>>();
    Ok(reaction.standard_enthalpy_change(&sources)?)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    let source_directory = match args.as_slice() {
        [] => None,
        [flag, directory] if flag == "--thermochemistry" => Some(Path::new(directory)),
        [flag] if flag == "--help" => {
            println!("{USAGE}");
            println!(
                "Compile the five USGS packs with xtask matdb-pack; see MINERAL_THERMOCHEMISTRY.md."
            );
            return Ok(());
        }
        _ => return Err(USAGE.into()),
    };
    let card = ConstitutiveModelCard {
        law: LawId(BILINEAR_BATCH_LAW.into()),
        law_version: 1,
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        parameters: BTreeMap::from([(
            RATE_PARAMETER.into(),
            LawParameter {
                value: 0.005,
                dims: BilinearRateConstant::DIMS,
            },
        )]),
        validity: ValidityDomain::unconstrained()
            .with_quantity(
                ELAPSED_TIME_AXIS,
                QuantitySpec::dimensional(Time::DIMS),
                0.0,
                100.0,
            )
            .with_quantity(
                VOLUME_AXIS,
                QuantitySpec::dimensional(Volume::DIMS),
                1.0,
                1.0,
            )
            .with_quantity(
                TEMPERATURE_AXIS,
                QuantitySpec::dimensional(Temperature::DIMS),
                298.15,
                298.15,
            )
            .with_quantity(
                REACTANT_A_CONCENTRATION_AXIS,
                QuantitySpec::dimensional(AmountConcentration::DIMS),
                0.0,
                10.0,
            )
            .with_quantity(
                REACTANT_B_CONCENTRATION_AXIS,
                QuantitySpec::dimensional(AmountConcentration::DIMS),
                0.0,
                20.0,
            ),
        sources: vec![hash_domain(
            "mineral-leaching-authored-example-v1",
            b"synthetic k=0.005 m3/(mol s); dolomite bicarbonate endpoint; V=1 m3; T=298.15 K",
        )],
        provenance: Provenance {
            source: "authored demonstrator; no measured mineral rate or calibration claim".into(),
            license: "project example".into(),
            artifact: None,
        },
    };
    let mineral = SpeciesId::new("CaMg(CO3)2(s)")?;
    let acid = SpeciesId::new("H+(aq)")?;
    let model = BilinearBatchReaction::from_card(
        MineralReaction::DolomiteBicarbonateNeutralization.build()?,
        mineral.clone(),
        acid.clone(),
        card,
    )?;
    println!("Authored rate-law demonstration; not calibrated mineral kinetics.");
    println!("Inventory/volume is empirical availability, not solid activity or surface area.");
    println!("Reaction: {}", model.reaction().id().as_str());
    let heat = source_directory
        .map(|directory| reference_enthalpy(directory, model.reaction()))
        .transpose()?;
    for initial_acid in [20.0, 2.0] {
        let initial = model
            .reaction()
            .terms()
            .iter()
            .map(|term| {
                let mol = if term.species == mineral {
                    10.0
                } else if term.species == acid {
                    initial_acid
                } else {
                    0.0
                };
                (term.species.clone(), Amount::new(mol))
            })
            .collect::<Vec<_>>();
        let result = model.advance(
            &initial,
            Volume::new(1.0),
            Time::new(10.0),
            Temperature::new(298.15),
            &QueryPoint::new(),
        )?;
        println!(
            "\nInitial H+: {initial_acid:.3} mol; contact time: {:.3} s",
            result.elapsed.value()
        );
        println!(
            "Actual extent: {:.9} mol equation; availability bound: {:.9} mol equation",
            result.extent.value(),
            result.limiting_extent.value()
        );
        for account in &result.reactants {
            println!(
                "  {}: consumed {:.9} mol; remaining {:.9} mol",
                account.species.as_str(),
                account.consumed.value(),
                account.remaining.value()
            );
        }
        println!("Complete final inventory:");
        for (species, amount) in &result.final_inventory {
            println!("  {}: {:.9} mol", species.as_str(), amount.value());
        }
        if let Some(heat) = &heat {
            println!(
                "Standard reaction enthalpy [J/mol extent]: {}",
                heat.molar().value()
            );
            println!(
                "Reaction enthalpy at computed extent [J]: {:.9}",
                heat.at_extent(result.extent)?.value()
            );
            println!(
                "298.15 K standard-state heat; mixing, solution activities and temperature evolution are separate."
            );
        }
        println!(
            "Original source: {}",
            result.model().source().provenance.source
        );
        println!(
            "Exact equation conservation: {:?}",
            result.model().reaction().conservation()
        );
    }
    Ok(())
}
