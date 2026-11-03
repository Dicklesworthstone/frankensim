//! Ideal-stage extraction followed by stripping using authored equilibrium fixtures.
//!
//! Run: `cargo run -p fs-material --example mining_countercurrent_sx`.
//! The numerical tables are explicitly synthetic, not copper/extractant data.
//! A process study must supply original measured isotherms and source contexts.

use fs_blake3::ContentHash;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance, QueryPoint,
};
use fs_material::solvent_extraction::{
    AQUEOUS_CONCENTRATION_AXIS, AQUEOUS_KNOT_PREFIX, ISOTHERM_LAW, ORGANIC_CONCENTRATION_AXIS,
    ORGANIC_KNOT_PREFIX, SolventExtractionIsotherm, TEMPERATURE_AXIS, VOLUME_RATIO_AXIS,
    countercurrent::{
        AmountFlowRate, CascadeNumerics, CountercurrentCascade, PhaseFlow, StageSpecification,
    },
};
use fs_qty::{
    AmountConcentration, Dimensionless, QuantitySpec, Temperature, VolumetricFlowRate,
    chemistry::ElementId,
};
use std::collections::BTreeMap;

fn fixture(
    slope: f64,
    maximum: f64,
    label: &str,
) -> Result<StageSpecification, Box<dyn std::error::Error>> {
    let mut parameters = BTreeMap::new();
    for (index, &(x, y)) in [(0.0, 0.0), (maximum, slope * maximum)].iter().enumerate() {
        for (prefix, value) in [(AQUEOUS_KNOT_PREFIX, x), (ORGANIC_KNOT_PREFIX, y)] {
            parameters.insert(
                format!("{prefix}{index}"),
                LawParameter {
                    value,
                    dims: AmountConcentration::DIMS,
                },
            );
        }
    }
    let source = ConstitutiveModelCard {
        law: LawId(ISOTHERM_LAW.into()),
        law_version: 1,
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        parameters,
        validity: ValidityDomain::unconstrained()
            .with_quantity(
                AQUEOUS_CONCENTRATION_AXIS,
                QuantitySpec::dimensional(AmountConcentration::DIMS),
                0.0,
                maximum,
            )
            .with_quantity(
                ORGANIC_CONCENTRATION_AXIS,
                QuantitySpec::dimensional(AmountConcentration::DIMS),
                0.0,
                slope * maximum,
            )
            .with_quantity(
                VOLUME_RATIO_AXIS,
                QuantitySpec::dimensional(Dimensionless::DIMS),
                0.1,
                10.0,
            )
            .with_quantity(
                TEMPERATURE_AXIS,
                QuantitySpec::dimensional(Temperature::DIMS),
                298.15,
                298.15,
            ),
        sources: vec![ContentHash([93; 32])],
        provenance: Provenance {
            source: format!("authored {label} balance fixture; no ore or extractant calibration"),
            license: "project example fixture".into(),
            artifact: Some(ContentHash([93; 32])),
        },
    };
    Ok(StageSpecification::new(
        SolventExtractionIsotherm::from_card(ElementId::new("Cu")?, source)?,
        Temperature::new(298.15),
        QueryPoint::new(),
    )?)
}
fn flow(volume: f64, concentration: f64) -> Result<PhaseFlow, Box<dyn std::error::Error>> {
    Ok(PhaseFlow::new(
        ElementId::new("Cu")?,
        VolumetricFlowRate::new(volume),
        AmountConcentration::new(concentration),
    )?)
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let numerics = CascadeNumerics {
        absolute_flow_tolerance: AmountFlowRate::new(1e-12),
        relative_flow_tolerance: Dimensionless::new(1e-12),
        max_root_evaluations: 160,
        max_stage_evaluations: 20_000,
    };
    let aqueous_feed = flow(2.0, 6.0)?;
    let barren_organic = flow(1.0, 0.0)?;
    let extraction_source = fixture(4.0, 16.0, "extraction")?;
    println!("Authored ideal-equilibrium fixture; no plant or extractant calibration.");
    println!("Fixed feeds: aqueous=2 m3/s at 6 mol/m3; organic=1 m3/s at 0 mol/m3.");
    println!("stages,raffinate_mol_m3,organic_out_mol_m3,transfer_mol_s,residual_mol_s");
    for count in 1..=4 {
        let train = CountercurrentCascade::new(vec![extraction_source.clone(); count])?;
        let result = train.solve(&aqueous_feed, &barren_organic, numerics, || false)?;
        println!(
            "{count},{:.9},{:.9},{:.9},{:.3e}",
            result.aqueous_outlet.concentration().value(),
            result.organic_outlet.concentration().value(),
            result.transfer_to_organic.value(),
            result.conservation_residual.value()
        );
    }
    let extraction = CountercurrentCascade::new(vec![extraction_source; 3])?.solve(
        &aqueous_feed,
        &barren_organic,
        numerics,
        || false,
    )?;
    // The loaded organic is reused as a real stream, with its actual flow and stock.
    let stripping = CountercurrentCascade::new(vec![fixture(0.25, 64.0, "stripping")?; 3])?.solve(
        &flow(0.25, 0.0)?,
        &extraction.organic_outlet,
        numerics,
        || false,
    )?;
    println!("\nThree-stage stripping of the three-stage extraction's organic outlet:");
    println!(
        "aqueous_product_mol_m3={:.9}; aqueous_product_mol_s={:.9}; organic_residual_mol_s={:.9}",
        stripping.aqueous_outlet.concentration().value(),
        stripping.aqueous_outlet.amount_flow().value(),
        stripping.organic_outlet.amount_flow().value()
    );
    println!(
        "stage,aqueous_out_mol_m3,organic_out_mol_m3,transfer_to_organic_mol_s,balance_residual_mol_s"
    );
    for (index, stage) in stripping.stages.iter().enumerate() {
        println!(
            "{},{:.9},{:.9},{:.9},{:.3e}",
            index + 1,
            stage.aqueous_outlet.concentration().value(),
            stage.organic_outlet.concentration().value(),
            stage.transfer_to_organic.value(),
            stage.conservation_residual.value()
        );
    }
    println!("Each stage retains the original source and its actual O/A/temperature context.");
    Ok(())
}
