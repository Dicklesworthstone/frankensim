//! Run an authored Bond grinding and additive-volume slurry scenario with real packs.
//!
//! Usage: `cargo run -p fs-material --example mineral_slurry -- \
//! magnetite.fsmatpk pyrite.fsmatpk dolomite.fsmatpk water.fsmatpk [SOLIDS_MASS_FRACTION]`
//!
//! The blend is 50% magnetite, 20% pyrite and 30% dolomite by dry mass.
//! Work index, power, sizes, utilization and initial water are authored inputs.
//! The optional target is a fraction, not percent; its default is 0.65.
//! The resulting volume flow also drives two separately authored continuous
//! recovery cards, with explicit active volume and a 3 g/t contained-metal assay.

use std::{collections::BTreeMap, error::Error, path::Path};

use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, MatDbError, NormalizedPack,
    Provenance, QueryPoint, SelectionPolicy,
};
use fs_material::{
    comminution::{BondGrindingCard, BondSizeDomain, BondTestContext, BondWorkIndex},
    recovery::{
        MAXIMUM_RATE_PARAMETER, RATE_PARAMETER, RecoveryError, ULTIMATE_RECOVERY_PARAMETER,
        continuous::{
            CONTINUOUS_SINGLE_RATE_LAW, CONTINUOUS_UNIFORM_RATE_LAW, ContinuousRecoveryCard,
            MEAN_RESIDENCE_TIME_AXIS, TANK_COUNT_AXIS,
        },
    },
    slurry::{DensitySource, SlurryComponent, SlurryStream, SolidsMassFraction},
};
use fs_qty::{
    Dimensionless, Dims, Frequency, Length, MassFlowRate, Power, QuantitySpec, Time, Volume,
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

fn authored_treatment(uniform: bool) -> Result<ContinuousRecoveryCard, RecoveryError> {
    let (law, parameter, rate) = if uniform {
        (
            CONTINUOUS_UNIFORM_RATE_LAW,
            MAXIMUM_RATE_PARAMETER,
            2.0 / 21600.0,
        )
    } else {
        (CONTINUOUS_SINGLE_RATE_LAW, RATE_PARAMETER, 1.0 / 21600.0)
    };
    let declaration = format!(
        "Authored stationary local law, no measured ore calibration: {law}, version 1, \
         Rinf=0.9; {parameter}={rate:.17e} s^-1; equal full ideal tanks; constant local \
         rates and carrier flow; no hidden state. Mean-residence support 43200..172800 s, \
         tank count 1..16, pH 9..11. Separate declaration from finite batch test evidence."
    );
    let artifact = hash_domain(
        "mineral-slurry-authored-continuous-v1",
        declaration.as_bytes(),
    );
    ContinuousRecoveryCard::from_card(ConstitutiveModelCard {
        law: LawId(law.into()),
        law_version: 1,
        parameters: BTreeMap::from([
            (
                ULTIMATE_RECOVERY_PARAMETER.into(),
                LawParameter {
                    value: 0.9,
                    dims: Dims::NONE,
                },
            ),
            (
                parameter.into(),
                LawParameter {
                    value: rate,
                    dims: Frequency::DIMS,
                },
            ),
        ]),
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained()
            .with_quantity(
                MEAN_RESIDENCE_TIME_AXIS,
                QuantitySpec::dimensional(Time::DIMS),
                43200.0,
                172800.0,
            )
            .with_quantity(
                TANK_COUNT_AXIS,
                QuantitySpec::dimensional(Dims::NONE),
                1.0,
                16.0,
            )
            .with_quantity("pH", QuantitySpec::dimensional(Dims::NONE), 9.0, 11.0),
        sources: vec![artifact],
        provenance: Provenance {
            source: declaration,
            license: "project-owned synthetic example".into(),
            artifact: Some(artifact),
        },
    })
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
        "Authored downstream treatment: 8 equal tanks, 17000 m3 active volume, pH=10, 3 g/t contained-metal feed assay."
    );
    println!(
        "Both separately authored continuous cards have Rinf=0.9 and mean local rate 1/(6 h); uniform rates span 0..2/(6 h)."
    );
    // The assay is an independent authored input, not inferred from a mineral
    // formula or density source. 3 g/dry metric tonne = 3e-6 kg/kg dry solids.
    let feed_metal = MassFlowRate::new(state.dry_solids_mass_flow.value() * 3.0e-6);
    let treatment_conditions =
        QueryPoint::new().with_quantity("pH", QuantitySpec::dimensional(Dims::NONE), 10.0)?;
    for (label, uniform) in [("single", false), ("uniform", true)] {
        let model = authored_treatment(uniform)?;
        match model.evaluate_train(
            8,
            Volume::new(17000.0),
            state.total_volume_flow,
            feed_metal,
            &treatment_conditions,
        ) {
            Ok(treatment) => println!(
                "Treatment {label}: actual_flow_m3_s={} mean_residence_h={} recovery={} recovered_metal_kg_h={} tail_metal_kg_h={} mass_residual_kg_s={} source_card={}",
                treatment.actual_train_flow.value(),
                treatment.mean_residence_time.value() / 3600.0,
                treatment.recovery.value(),
                treatment.recovered_metal_flow.value() * 3600.0,
                treatment.tail_metal_flow.value() * 3600.0,
                treatment.conservation_residual.value(),
                treatment.source_card.content_hash(),
            ),
            Err(RecoveryError::OutsideValidity) => println!(
                "Treatment {label}: OutsideValidity; this slurry flow gives a mean residence outside the independently declared 12..48 h support."
            ),
            Err(error) => return Err(error.into()),
        }
    }
    println!(
        "Additive constituent volumes exclude gas, contraction, dissolution and settling. Source uncertainty is retained without a propagated slurry-density error bound."
    );
    Ok(())
}
