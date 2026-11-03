//! Run: cargo run -p fs-material --example mineral_solvent_extraction
//! Optional arguments: FEED_CU_MOL PASSED_CHARGE_C (defaults: 100, 3600000).
//!
//! Authored synthetic extraction → stripping → productive copper EW. These
//! isotherms are illustrative input data, not measurements or LIX calibration.
//! Both source tables use organic C versus aqueous C, in mol/m³. The operator
//! accounts for analytical Cu only; the example explicitly assigns strip Cu
//! to Cu(II) and supplies water/protons before invoking the separate EW owner.

use std::{collections::BTreeMap, error::Error};

use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance, QueryPoint,
};
use fs_material::solvent_extraction::{
    AQUEOUS_CONCENTRATION_AXIS, AQUEOUS_KNOT_PREFIX, ISOTHERM_LAW, ORGANIC_CONCENTRATION_AXIS,
    ORGANIC_KNOT_PREFIX, PhaseInventory, SolventExtractionIsotherm, TEMPERATURE_AXIS,
    VOLUME_RATIO_AXIS,
};
use fs_qty::{
    Amount, AmountConcentration, Dimensionless, ElectricCharge, MolarMass, QuantitySpec,
    Temperature, Volume, chemistry::ElementId,
};
use fs_thermochem::{
    SpeciesId,
    stoichiometry::{electrolysis::ElectrolysisReaction, minerals::MineralReaction},
};

fn isotherm(
    name: &str,
    chemistry: f64,
    knots: &[(f64, f64)],
) -> Result<(SolventExtractionIsotherm, QueryPoint), Box<dyn Error>> {
    let description = format!(
        "Authored synthetic {name} equilibrium table: {knots:?} mol/m3 (aqueous, organic); \
         Cu analytical basis; fixed 298.15 K; O/A 0.1..10; chemistry scenario {chemistry}. \
         Not measured extractant performance. No acid, extractant or phase-volume balance."
    );
    let artifact = hash_domain(
        "mineral-solvent-extraction-authored-source-v1",
        description.as_bytes(),
    );
    let mut parameters = BTreeMap::new();
    for (i, &(aqueous, organic)) in knots.iter().enumerate() {
        for (prefix, value) in [
            (AQUEOUS_KNOT_PREFIX, aqueous),
            (ORGANIC_KNOT_PREFIX, organic),
        ] {
            parameters.insert(
                format!("{prefix}{i}"),
                LawParameter {
                    value,
                    dims: AmountConcentration::DIMS,
                },
            );
        }
    }
    let q = QuantitySpec::dimensional;
    let source = ConstitutiveModelCard {
        law: LawId(ISOTHERM_LAW.into()),
        law_version: 1,
        parameters,
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained()
            .with_quantity(
                AQUEOUS_CONCENTRATION_AXIS,
                q(AmountConcentration::DIMS),
                knots[0].0,
                knots.last().unwrap().0,
            )
            .with_quantity(
                ORGANIC_CONCENTRATION_AXIS,
                q(AmountConcentration::DIMS),
                knots[0].1,
                knots.last().unwrap().1,
            )
            .with_quantity(VOLUME_RATIO_AXIS, q(Dimensionless::DIMS), 0.1, 10.0)
            .with_quantity(TEMPERATURE_AXIS, q(Temperature::DIMS), 298.15, 298.15)
            .with_quantity(
                "authored-chemistry",
                q(Dimensionless::DIMS),
                chemistry,
                chemistry,
            ),
        sources: vec![artifact],
        provenance: Provenance {
            source: description,
            license: "project-authored example".into(),
            artifact: Some(artifact),
        },
    };
    let conditions =
        QueryPoint::new().with_quantity("authored-chemistry", q(Dimensionless::DIMS), chemistry)?;
    Ok((
        SolventExtractionIsotherm::from_card(ElementId::new("Cu")?, source)?,
        conditions,
    ))
}

fn phase(volume: f64, amount: f64) -> Result<PhaseInventory, Box<dyn Error>> {
    Ok(PhaseInventory::new(
        ElementId::new("Cu")?,
        Volume::new(volume),
        Amount::new(amount),
    )?)
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.len() == 1 && (args[0] == "--help" || args[0] == "-h") {
        println!(
            "Usage: mineral_solvent_extraction [FEED_CU_MOL PASSED_CHARGE_C]\n\
            Defaults: 100 mol feed Cu, 3600000 C, productive current efficiency 0.95.\n\
            1 m3 aqueous feed contacts 0.5 m3 unloaded organic; the loaded organic\n\
            then contacts 0.1 m3 Cu-free strip solution using a separate isotherm.\n\
            Synthetic source tables refuse unsupported equilibria; EW refuses stock overdraw.\n\
            No commercial extractant, acid balance, actual stage efficiency or cell voltage is inferred."
        );
        return Ok(());
    }
    if !args.is_empty() && args.len() != 2 {
        return Err("expected FEED_CU_MOL PASSED_CHARGE_C, or no arguments; use --help".into());
    }
    let feed = args.first().map_or(Ok(100.0), |s| s.parse::<f64>())?;
    let charge = args.get(1).map_or(Ok(3_600_000.0), |s| s.parse::<f64>())?;
    let (extraction, extraction_conditions) = isotherm(
        "extraction",
        1.0,
        &[(0.0, 0.0), (50.0, 200.0), (100.0, 300.0)],
    )?;
    let (stripping, stripping_conditions) = isotherm(
        "stripping",
        2.0,
        &[(0.0, 0.0), (400.0, 20.0), (800.0, 40.0)],
    )?;
    let extracted = extraction.equilibrate(
        &phase(1.0, feed)?,
        &phase(0.5, 0.0)?,
        Temperature::new(298.15),
        &extraction_conditions,
    )?;
    let stripped = stripping.equilibrate(
        &phase(0.1, 0.0)?,
        &extracted.organic,
        Temperature::new(298.15),
        &stripping_conditions,
    )?;

    // This authored species assignment is explicit: SX's analytical Cu balance
    // does not itself establish Cu(II), free acid or solvent water inventory.
    let inventory = [
        ("Cu[2+](aq)", stripped.aqueous.amount().value()),
        ("H2O(l)", 5_550.0),
        ("H+(aq)", 100.0),
        ("Cu(s)", 0.0),
        ("O2(g)", 0.0),
    ]
    .into_iter()
    .map(|(name, n)| Ok((SpeciesId::new(name)?, Amount::new(n))))
    .collect::<Result<Vec<_>, Box<dyn Error>>>()?;
    let ew = ElectrolysisReaction::new(
        MineralReaction::CopperElectrowinning.build()?,
        SpeciesId::new("Cu(s)")?,
        MolarMass::new(0.063546),
        2,
    )?
    .advance_charge(
        &inventory,
        ElectricCharge::new(charge),
        Dimensionless::new(0.95),
    )?;
    let remaining_copper = ew
        .modeled_final_inventory
        .iter()
        .find(|(id, _)| id.as_str() == "Cu[2+](aq)")
        .ok_or("missing electrolyte Cu")?
        .1
        .value();
    println!(
        "Authored extraction/stripping tables; ideal equilibrium contacts, not a calibrated plant."
    );
    println!(
        "Extraction source: {}",
        extraction.source().sources[0].to_hex()
    );
    println!(
        "Stripping source: {}",
        stripping.source().sources[0].to_hex()
    );
    println!("Feed Cu [mol]: {feed}");
    println!(
        "Extraction raffinate Cu [mol]: {}",
        extracted.aqueous.amount().value()
    );
    println!(
        "Loaded organic Cu [mol]: {}",
        extracted.organic.amount().value()
    );
    println!(
        "Strip aqueous Cu [mol]: {}",
        stripped.aqueous.amount().value()
    );
    println!(
        "Strip aqueous Cu [mol/m3]: {}",
        stripped.aqueous.concentration().value()
    );
    println!(
        "Stripped organic Cu [mol]: {}",
        stripped.organic.amount().value()
    );
    println!(
        "Extraction Cu residual [mol]: {}",
        extracted.conservation_residual.value()
    );
    println!(
        "Stripping Cu residual [mol]: {}",
        stripped.conservation_residual.value()
    );
    println!(
        "EW Cu deposited [mol]: {}",
        ew.product_yield.amount().value()
    );
    println!("EW Cu deposited [kg]: {}", ew.product_yield.mass().value());
    println!("EW electrolyte Cu remaining [mol]: {remaining_copper}");
    println!("Unallocated charge [C]: {}", ew.unallocated_charge.value());
    println!(
        "Overall Cu residual [mol]: {}",
        extracted.aqueous.amount().value()
            + stripped.organic.amount().value()
            + remaining_copper
            + ew.product_yield.amount().value()
            - feed
    );
    println!(
        "SX accounts for analytical Cu only. EW water/protons and Cu(II) assignment are authored inputs; unallocated-charge side chemistry remains unmodeled."
    );
    Ok(())
}
