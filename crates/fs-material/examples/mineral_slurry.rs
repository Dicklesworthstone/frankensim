//! Run an authored Bond grinding and additive-volume slurry scenario with real packs.
//!
//! Usage: `cargo run -p fs-material --example mineral_slurry -- \
//! magnetite.fsmatpk pyrite.fsmatpk dolomite.fsmatpk water.fsmatpk [SOLIDS_MASS_FRACTION]`
//!
//! The blend is 50% magnetite, 20% pyrite and 30% dolomite by dry mass.
//! Work index, power, sizes, utilization and initial water are authored inputs.
//! The optional target is a fraction, not percent; its default is 0.65.

use std::{error::Error, path::Path};

use fs_matdb::{MatDbError, NormalizedPack, Provenance, QueryPoint, SelectionPolicy};
use fs_material::{
    comminution::{BondGrindingCard, BondSizeDomain, BondTestContext, BondWorkIndex},
    slurry::{DensitySource, SlurryComponent, SlurryStream, SolidsMassFraction},
};
use fs_qty::{
    Dimensionless, Dims, Length, MassFlowRate, Power, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

const USAGE: &str = "usage: mineral_slurry MAGNETITE_PACK PYRITE_PACK DOLOMITE_PACK WATER_PACK [SOLIDS_MASS_FRACTION, default 0.65]";

fn load(path: &Path, expected_id: &str) -> Result<NormalizedPack, Box<dyn Error>> {
    let bytes = std::fs::read(path)?;
    let decoded = NormalizedPack::from_bytes(&bytes)?;
    let pack = NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes)?;
    if pack.to_bytes() != bytes {
        return Err(format!("{} is not a canonical material pack", path.display()).into());
    }
    if pack.pack_id() != expected_id {
        return Err(format!("expected pack {expected_id}, received {}", pack.pack_id()).into());
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

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args_os().skip(1).collect();
    if args.len() == 1 && args[0] == "--help" {
        println!("{USAGE}");
        return Ok(());
    }
    if !(4..=5).contains(&args.len()) {
        return Err(USAGE.into());
    }
    let target = match args.get(4) {
        Some(value) => value
            .to_str()
            .ok_or("solids mass fraction must be a numeric string")?
            .parse::<f64>()?,
        None => 0.65,
    };
    let target = SolidsMassFraction::new(Dimensionless::new(target))?;
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

    let feed_size = Length::new(0.0025);
    let product_size = Length::new(0.0001);
    let grinding = BondGrindingCard::new(
        BondWorkIndex::from_kwh_per_metric_tonne(15.0)?,
        BondTestContext {
            specimen: "authored mineral blend; not a measured ore".into(),
            method: "authored Bond ball-mill scenario with prepared feed".into(),
            closing_screen: Length::new(150e-6),
            provenance: Provenance {
                source: "authored example inputs, not a site calibration".into(),
                license: "project example".into(),
                artifact: None,
            },
        },
        BondSizeDomain::new([feed_size; 2], [product_size; 2])?,
    )?;
    let energy = grinding.specific_energy(feed_size, product_size)?;
    let dry = grinding.dry_throughput(
        feed_size,
        product_size,
        Power::new(12e6),
        Dimensionless::new(0.9),
    )?;
    let reference = mineral_point()?;
    let mut solids = Vec::with_capacity(3);
    for ((id, fraction), pack) in [("magnetite", 0.5), ("pyrite", 0.2), ("dolomite", 0.3)]
        .into_iter()
        .zip(&packs[..3])
    {
        let density = DensitySource::resolve(
            pack,
            "density",
            &reference,
            SelectionPolicy::SingleClaimOnly,
        )?;
        solids.push(SlurryComponent::new(
            id,
            MassFlowRate::new(dry.value() * fraction),
            density,
        )?);
    }
    let water = DensitySource::resolve(
        &packs[3],
        "density",
        &water_point()?,
        SelectionPolicy::SingleClaimOnly,
    )?;
    let initial = SlurryStream::new_additive_volume(
        solids,
        SlurryComponent::new("process-water", MassFlowRate::new(50.0), water)?,
    )?;
    let dilution = initial.dilute_to_mass_fraction(target)?;
    let state = dilution.stream.state();
    let density = state.density.ok_or("zero-flow slurry has no density")?;
    let mass_fraction = state
        .solids_mass_fraction
        .ok_or("zero-flow slurry has no solids mass fraction")?;
    let volume_fraction = state
        .solids_volume_fraction
        .ok_or("zero-flow slurry has no solids volume fraction")?;

    println!("Authored mineral-processing scenario; not a plant forecast.");
    println!("Dry blend: 50% magnetite, 20% pyrite, 30% dolomite by mass.");
    println!(
        "Inputs: Wi=15 kWh/dry metric tonne, F80=2500 um, P80=100 um, net pinion power=12 MW, utilization=0.9."
    );
    println!("Grinding energy [J/kg dry solids]: {}", energy.value());
    println!("Dry-solids throughput [kg/s]: {}", dry.value());
    println!(
        "Initial water [kg/s]: {}",
        initial.state().liquid_mass_flow.value()
    );
    println!("Added water [kg/s]: {}", dilution.added_liquid_flow.value());
    println!("Total water [kg/s]: {}", state.liquid_mass_flow.value());
    println!(
        "Total slurry mass flow [kg/s]: {}",
        state.total_mass_flow.value()
    );
    println!(
        "Total additive volume flow [m3/s]: {}",
        state.total_volume_flow.value()
    );
    println!("Pulp density [kg/m3]: {}", density.value());
    println!("Solids mass fraction [1]: {}", mass_fraction.value());
    println!("Solids volume fraction [1]: {}", volume_fraction.value());
    for component in dilution
        .stream
        .solids()
        .iter()
        .chain([dilution.stream.liquid()])
    {
        let source = component.density_source();
        println!(
            "Source {}: {} ({})",
            component.id(),
            source.pack_id(),
            source.pack_hash()
        );
    }
    println!(
        "Water density is queried at 298.15 K and 100000 Pa; mineral measurement conditions remain unstated."
    );
    println!(
        "Additive constituent volumes exclude gas, contraction, dissolution and settling. Source uncertainty is retained without a propagated slurry-density error bound."
    );
    Ok(())
}
