//! Closed rougher/cleaner recycle using the original mineral and water packs.
//!
//! `cargo run -p fs-material --example mineral_recycle -- \
//! MAGNETITE_PACK PYRITE_PACK DOLOMITE_PACK WATER_PACK [RECYCLE_FRACTION]`
//!
//! Partitions, feed rates and the default 0.8 cleaner-tail return are authored
//! examples. No measured flotation response, liberation or capacity is inferred.

use std::{error::Error, path::Path};

use fs_blake3::hash_domain;
use fs_matdb::{MatDbError, NormalizedPack, Provenance, QueryPoint, SelectionPolicy};
use fs_material::slurry::{
    DensitySource, SlurryComponent, SlurryStream,
    recycle::RougherCleanerCircuit,
    separation::{ComponentPartition, SelectiveSeparator},
};
use fs_qty::{
    Dimensionless, Dims, MassFlowRate, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

const USAGE: &str = "usage: mineral_recycle MAGNETITE_PACK PYRITE_PACK DOLOMITE_PACK WATER_PACK [RECYCLE_FRACTION, default 0.8]";
const SOLIDS: [&str; 3] = ["magnetite", "pyrite", "dolomite"];

fn load(path: &Path, expected_id: &str) -> Result<NormalizedPack, Box<dyn Error>> {
    let bytes = std::fs::read(path)?;
    let decoded = NormalizedPack::from_bytes(&bytes)?;
    let pack = NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes)?;
    if pack.to_bytes() != bytes {
        return Err(format!("{} is not a canonical pack", path.display()).into());
    }
    if pack.pack_id() != expected_id {
        return Err(format!("expected {expected_id}, received {}", pack.pack_id()).into());
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

fn authored_source(text: String) -> Provenance {
    let artifact = hash_domain("mineral-recycle-authored-scenario-v1", text.as_bytes());
    Provenance {
        source: text,
        license: "project-owned synthetic example".into(),
        artifact: Some(artifact),
    }
}

fn separator(
    stage: &str,
    fractions: [f64; 3],
    water: f64,
) -> Result<SelectiveSeparator, Box<dyn Error>> {
    let source = authored_source(format!(
        "Authored {stage} separation, not measured plant performance: \
        component mass fractions to concentrate magnetite={:.17}, pyrite={:.17}, dolomite={:.17}; \
        independent water fraction={water:.17}. Fixed partitions assumed applicable to the recycled blend.",
        fractions[0], fractions[1], fractions[2]
    ));
    Ok(SelectiveSeparator::new(
        SOLIDS
            .into_iter()
            .zip(fractions)
            .map(|(id, fraction)| ComponentPartition::new(id, Dimensionless::new(fraction)))
            .collect::<Result<Vec<_>, _>>()?,
        Dimensionless::new(water),
        source,
    )?)
}

fn optional(value: Option<Dimensionless>) -> String {
    value
        .map(|value| format!("{:.12}", value.value()))
        .unwrap_or_default()
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args_os().skip(1).collect();
    if args.len() == 1 && args[0] == "--help" {
        println!("{USAGE}");
        return Ok(());
    }
    if !(4..=5).contains(&args.len()) {
        return Err(USAGE.into());
    }
    let returned = match args.get(4) {
        Some(value) => value
            .to_str()
            .ok_or("return fraction must be a numeric string")?
            .parse::<f64>()?,
        None => 0.8,
    };
    let packs = [
        load(
            Path::new(&args[0]),
            "magnetite-handbook-mineralogy-reference",
        )?,
        load(Path::new(&args[1]), "pyrite-handbook-mineralogy-reference")?,
        load(
            Path::new(&args[2]),
            "dolomite-handbook-mineralogy-reference",
        )?,
        load(Path::new(&args[3]), "water-liquid-iapws-sr6-08-0.1mpa")?,
    ];
    let point = mineral_point()?;
    let mut solids = Vec::with_capacity(3);
    for ((id, flow), pack) in SOLIDS.into_iter().zip([100.0, 20.0, 80.0]).zip(&packs[..3]) {
        solids.push(SlurryComponent::new(
            id,
            MassFlowRate::new(flow),
            DensitySource::resolve(pack, "density", &point, SelectionPolicy::SingleClaimOnly)?,
        )?);
    }
    let water = SlurryComponent::new(
        "water",
        MassFlowRate::new(100.0),
        DensitySource::resolve(
            &packs[3],
            "density",
            &water_point()?,
            SelectionPolicy::SingleClaimOnly,
        )?,
    )?;
    let fresh = SlurryStream::new_additive_volume(solids, water)?;
    let circuit = RougherCleanerCircuit::new(
        separator("rougher", [0.95, 0.3, 0.05], 0.2)?,
        separator("cleaner", [0.9, 0.2, 0.1], 0.25)?,
        Dimensionless::new(returned),
        authored_source(format!(
            "Authored closed rougher-cleaner circuit: fresh mass flows 100/20/80 kg/s \
            magnetite/pyrite/dolomite and 100 kg/s water; nonselective fraction {returned:.17} of cleaner \
            tails returns to rougher feed, complement purged. No size change, reaction or extra makeup."
        )),
    )?;
    let result = circuit.solve(&fresh)?;
    println!(
        "# Authored steady rougher/cleaner circuit with cleaner-tail return={returned}; remaining cleaner tail is purged."
    );
    println!(
        "# Prescribed stage partitions are assumed constant for the recycled blend; no calibration or capacity qualification is asserted."
    );
    println!(
        "# All density claims and query contexts are retained. Mineral reference measurement conditions remain unstated; water is 298.15 K and 100000 Pa."
    );
    for component in fresh.solids().iter().chain([fresh.liquid()]) {
        println!(
            "# Density source {}: {} ({})",
            component.id(),
            component.density_source().pack_id(),
            component.density_source().pack_hash()
        );
    }
    println!("stream,dry_solids_kg_s,water_kg_s,total_mass_kg_s,volume_m3_s,density_kg_m3");
    for (name, stream) in [
        ("fresh", &result.fresh_feed),
        ("mixed_rougher_feed", &result.mixed_feed),
        ("rougher_concentrate", &result.rougher.concentrate),
        ("rougher_tailings", &result.rougher.tailings),
        ("final_product", &result.cleaner.concentrate),
        ("cleaner_tailings", &result.cleaner.tailings),
        ("recycle", &result.recycle),
        ("purge", &result.purge),
        ("combined_external_tailings", &result.combined_tailings),
    ] {
        let state = stream.state();
        let density = state
            .density
            .map(|v| format!("{:.9}", v.value()))
            .unwrap_or_default();
        println!(
            "{name},{:.12},{:.12},{:.12},{:.12},{density}",
            state.dry_solids_mass_flow.value(),
            state.liquid_mass_flow.value(),
            state.total_mass_flow.value(),
            state.total_volume_flow.value()
        );
    }
    println!(
        "# Fractions below concern each named component mass, not elemental assays. Water is reported separately from dry solids."
    );
    println!(
        "component,fresh_kg_s,mixed_kg_s,product_kg_s,recycle_kg_s,purge_kg_s,circulating_load,overall_recovery,external_residual_kg_s,recycle_residual_kg_s"
    );
    for balance in result.components.iter().chain([&result.carrier]) {
        println!(
            "{},{:.12},{:.12},{:.12},{:.12},{:.12},{},{},{:.6e},{:.6e}",
            balance.component,
            balance.fresh_mass_flow.value(),
            balance.mixed_mass_flow.value(),
            balance.product_mass_flow.value(),
            balance.recycle_mass_flow.value(),
            balance.purge_mass_flow.value(),
            optional(balance.circulating_load),
            optional(balance.overall_recovery),
            balance.external_conservation_residual.value(),
            balance.recycle_closure_residual.value()
        );
    }
    println!("final_product_component,dry_mass_fraction");
    for component in &result.cleaner.components {
        println!(
            "{},{}",
            component.balance.component,
            optional(component.concentrate_dry_fraction)
        );
    }
    Ok(())
}
