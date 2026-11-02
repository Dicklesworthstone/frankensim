//! Prescribed rougher/cleaner separation with real density packs and closed outlets.
//!
//! Usage: `cargo run -p fs-material --example mineral_partition_circuit -- \
//! magnetite.fsmatpk pyrite.fsmatpk dolomite.fsmatpk water.fsmatpk`
//!
//! Feed is 100 kg/s magnetite, 20 pyrite, 80 dolomite and 100 water. Partitions
//! are authored scenarios, not measured flotation or cyclone performance. Dry
//! component fractions refer to these named solids, not elemental Fe/Cu assays.

use std::{error::Error, path::Path};

use fs_matdb::{MatDbError, NormalizedPack, Provenance, QueryPoint, SelectionPolicy};
use fs_material::slurry::{
    DensitySource, SlurryComponent, SlurryStream,
    separation::{ComponentPartition, SelectiveSeparator, SeparationError},
};
use fs_qty::{
    Dimensionless, Dims, MassFlowRate, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

const USAGE: &str =
    "usage: mineral_partition_circuit MAGNETITE_PACK PYRITE_PACK DOLOMITE_PACK WATER_PACK";
const SOLIDS: [&str; 3] = ["magnetite", "pyrite", "dolomite"];

fn load(path: &Path, expected_id: &str) -> Result<NormalizedPack, Box<dyn Error>> {
    let bytes = std::fs::read(path)?;
    let decoded = NormalizedPack::from_bytes(&bytes)?;
    let pack = NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes)?;
    if pack.pack_id() != expected_id || pack.to_bytes() != bytes {
        return Err(format!(
            "{} must be the canonical {expected_id} pack",
            path.display()
        )
        .into());
    }
    Ok(pack)
}

fn mineral_point() -> Result<QueryPoint, MatDbError> {
    let mut point = QueryPoint::new();
    for (axis, value) in [
        ("mineral-reference-model", 1.0),
        ("source-test-temperature-known", 0.0),
        ("source-test-pressure-known", 0.0),
        ("bulk-ore-qualified", 0.0),
    ] {
        point = point.with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)?;
    }
    Ok(point)
}

fn water_point() -> Result<QueryPoint, MatDbError> {
    QueryPoint::new()
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
        )?
        .with_quantity("phase-liquid", QuantitySpec::dimensional(Dims::NONE), 1.0)
}

fn stage(
    name: &str,
    solid_partitions: [f64; 3],
    water: f64,
) -> Result<SelectiveSeparator, SeparationError> {
    SelectiveSeparator::new(
        SOLIDS
            .into_iter()
            .zip(solid_partitions)
            .map(|(id, fraction)| ComponentPartition::new(id, Dimensionless::new(fraction)))
            .collect::<Result<Vec<_>, _>>()?,
        Dimensionless::new(water),
        Provenance {
            source: format!(
                "authored {name} partitions; no site calibration or performance measurement"
            ),
            license: "project example".into(),
            artifact: None,
        },
    )
}

fn report(label: &str, stream: &SlurryStream) {
    let state = stream.state();
    println!(
        "{label}: dry_kg_s={:.6} water_kg_s={:.6} volume_m3_s={:.9} density_kg_m3={:?}",
        state.dry_solids_mass_flow.value(),
        state.liquid_mass_flow.value(),
        state.total_volume_flow.value(),
        state.density.map(|value| value.value()),
    );
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args_os().skip(1).collect();
    if args.len() == 1 && args[0] == "--help" {
        println!("{USAGE}");
        return Ok(());
    }
    if args.len() != 4 {
        return Err(USAGE.into());
    }
    let mut solids = Vec::with_capacity(3);
    let mineral_query = mineral_point()?;
    for ((id, mass), path) in SOLIDS.into_iter().zip([100.0, 20.0, 80.0]).zip(&args[..3]) {
        let pack = load(
            Path::new(path),
            &format!("{id}-handbook-mineralogy-reference"),
        )?;
        let source = DensitySource::resolve(
            &pack,
            "density",
            &mineral_query,
            SelectionPolicy::SingleClaimOnly,
        )?;
        println!(
            "source {id}: pack={} density_kg_m3={} receipt={}",
            source.pack_id(),
            source.density().value(),
            source.answer().receipt.content_hash().to_hex()
        );
        solids.push(SlurryComponent::new(id, MassFlowRate::new(mass), source)?);
    }
    let water = load(Path::new(&args[3]), "water-liquid-iapws-sr6-08-0.1mpa")?;
    let source = DensitySource::resolve(
        &water,
        "density",
        &water_point()?,
        SelectionPolicy::SingleClaimOnly,
    )?;
    println!(
        "source water: pack={} density_kg_m3={} receipt={}",
        source.pack_id(),
        source.density().value(),
        source.answer().receipt.content_hash().to_hex()
    );
    let feed = SlurryStream::new_additive_volume(
        solids,
        SlurryComponent::new("water", MassFlowRate::new(100.0), source)?,
    )?;
    let rougher = stage("rougher", [0.95, 0.30, 0.05], 0.20)?.separate(&feed)?;
    let cleaner = stage("cleaner", [0.90, 0.20, 0.10], 0.25)?.separate(&rougher.concentrate)?;
    let combined_tailings = SlurryStream::mix(&[&rougher.tailings, &cleaner.tailings])?;
    report("feed", &feed);
    report("rougher concentrate", &rougher.concentrate);
    report("final concentrate", &cleaner.concentrate);
    report("combined tailings", &combined_tailings);
    for (original, (concentrate, tailings)) in feed.solids().iter().zip(
        cleaner
            .concentrate
            .solids()
            .iter()
            .zip(combined_tailings.solids()),
    ) {
        let recovery = concentrate.mass_flow().value() / original.mass_flow().value();
        let dry_fraction = concentrate.mass_flow().value()
            / cleaner.concentrate.state().dry_solids_mass_flow.value();
        let residual = concentrate.mass_flow().value() + tailings.mass_flow().value()
            - original.mass_flow().value();
        println!(
            "{}: final_dry_fraction={dry_fraction:.9} overall_recovery={recovery:.9} mass_residual_kg_s={residual:.3e}",
            original.id()
        );
    }
    let water_residual = cleaner.concentrate.liquid().mass_flow().value()
        + combined_tailings.liquid().mass_flow().value()
        - feed.liquid().mass_flow().value();
    println!("water_mass_residual_kg_s={water_residual:.3e}");
    println!(
        "Authored partitions; additive constituent volumes without gas. Reference grain densities do not qualify bulk ore or establish common temperature. No calibrated cyclone/flotation, liberation, or elemental assay is inferred."
    );
    Ok(())
}
