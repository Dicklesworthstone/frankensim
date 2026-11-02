//! Authored grinding, classification and component-specific recovery with balances.
//!
//! Run: `cargo run -p fs-material --example mineral_separation -- [RESIDENCE_SECONDS]`
//! The default is 120 seconds; all parameters are illustrative, not ore testwork.
//! Underflow is reported separately. This single pass does not solve mill recycle.

use std::{collections::BTreeMap, error::Error};

use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance, QueryPoint,
};
use fs_material::{
    comminution::{BondGrindingCard, BondSizeDomain, BondTestContext, BondWorkIndex},
    recovery::distributed::{
        RESIDENCE_TIME_AXIS, RESIDENCE_TIME_QUANTITY, RateDistribution, RatePopulation,
        RecoveryCard, RecoveryReactor, RecoveryTestContext,
    },
    separation::{
        CLASSIFIER_LAW_ID, ClassifierCard, PARTICLE_DIAMETER_AXIS, ParticleClass, ParticleStream,
        Partition, SOLIDS_FRACTION_AXIS,
    },
};
use fs_qty::{
    Dimensionless, Dims, Frequency, Length, MassFlowRate, Power, QuantitySpec, Time,
    chemistry::ElementId,
};

const USAGE: &str = "usage: mineral_separation [RESIDENCE_SECONDS, default 120; supported 0..3600]";
const CLASSIFIER_SOURCE: &str = "Authored classifier: corrected d50=0.0001 m, sharpness=2, \
water recovery to underflow=0.2; particle sizes=0.000001..0.003 m; \
solids mass fraction=0.2..0.8. These are illustrative parameters, not measured data.";

fn provenance(source: &str) -> Provenance {
    Provenance {
        source: source.into(),
        license: "project authored example".into(),
        artifact: Some(hash_domain(
            "mineral-separation-example-v1",
            source.as_bytes(),
        )),
    }
}

fn classifier() -> Result<ClassifierCard, Box<dyn Error>> {
    let source = provenance(CLASSIFIER_SOURCE);
    let card = ConstitutiveModelCard {
        law: LawId(CLASSIFIER_LAW_ID.into()),
        law_version: 1,
        parameters: [
            ("corrected-cut-size", 100e-6, Length::DIMS),
            ("sharpness", 2.0, Dims::NONE),
            ("water-recovery-to-underflow", 0.2, Dims::NONE),
        ]
        .into_iter()
        .map(|(name, value, dims)| (name.into(), LawParameter { value, dims }))
        .collect(),
        state_schema_version: 1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained()
            .with_quantity(
                PARTICLE_DIAMETER_AXIS,
                QuantitySpec::dimensional(Length::DIMS),
                1e-6,
                0.003,
            )
            .with_quantity(
                SOLIDS_FRACTION_AXIS,
                QuantitySpec::dimensional(Dims::NONE),
                0.2,
                0.8,
            ),
        sources: vec![source.artifact.ok_or("missing authored source hash")?],
        provenance: source,
    };
    Ok(ClassifierCard::from_model_card(&card)?)
}

fn kinetics(valuable: bool) -> Result<RecoveryCard, Box<dyn Error>> {
    let rates = if valuable {
        RateDistribution::Discrete {
            populations: vec![
                RatePopulation {
                    weight: Dimensionless::new(0.6),
                    rate: Frequency::new(0.02),
                },
                RatePopulation {
                    weight: Dimensionless::new(0.4),
                    rate: Frequency::new(0.002),
                },
            ],
        }
    } else {
        RateDistribution::Single {
            rate: Frequency::new(0.005),
        }
    };
    let source = if valuable {
        "Authored valuable-class fit: Rinf=.95; relative populations(.6,.02/s),(.4,.002/s); two equal steady mixed tanks; total residence0..3600s; pH9..11. No measured calibration."
    } else {
        "Authored gangue-class fit: Rinf=.05; k=.005/s; two equal steady mixed tanks; total residence0..3600s; pH9..11. No measured calibration."
    };
    Ok(RecoveryCard::new(
        Dimensionless::new(if valuable { 0.95 } else { 0.05 }),
        rates,
        RecoveryReactor::EqualCstrCascade { tanks: 2 },
        RecoveryTestContext {
            specimen: if valuable {
                "authored valuable particle classes"
            } else {
                "authored gangue particle classes"
            }
            .into(),
            method: "illustrative constant-rate steady two-tank closure".into(),
            provenance: provenance(source),
        },
        ValidityDomain::unconstrained()
            .with_quantity(RESIDENCE_TIME_AXIS, RESIDENCE_TIME_QUANTITY, 0.0, 3600.0)
            .with_quantity("pH", QuantitySpec::dimensional(Dims::NONE), 9.0, 11.0),
    )?)
}

fn report(name: &str, stream: &ParticleStream, copper: &ElementId) -> Result<(), Box<dyn Error>> {
    let grade = stream.element_grade(copper)?.map_or_else(
        || "undefined (zero solids)".into(),
        |value| format!("{:.6}%", 100.0 * value.value()),
    );
    println!(
        "{name}: dry solids={:.9} kg/s; Cu={:.9} kg/s; water={:.9} kg/s; Cu grade={grade}",
        stream.dry_solids()?.value(),
        stream.element_flow(copper)?.value(),
        stream.water().value()
    );
    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.len() == 1 && args[0] == "--help" {
        println!("{USAGE}");
        return Ok(());
    }
    if args.len() > 1 {
        return Err(USAGE.into());
    }
    let time = Time::new(
        args.first()
            .map(|s| s.parse::<f64>())
            .transpose()?
            .unwrap_or(120.0),
    );
    let feed_size = Length::new(0.0025);
    // The explicitly authored two-bin product has its discrete P80 at 200 um.
    let product_size = Length::new(200e-6);
    let grinding = BondGrindingCard::new(
        BondWorkIndex::from_kwh_per_metric_tonne(15.0)?,
        BondTestContext {
            specimen: "authored feed; no measured work index".into(),
            method: "illustrative Bond ball-mill energy closure".into(),
            closing_screen: Length::new(250e-6),
            provenance: provenance(
                "Authored Bond inputs: Wi15 kWh/dry metric tonne; F80=2500um; P80=200um; net pinion power1MW; utilization.9.",
            ),
        },
        BondSizeDomain::new([feed_size; 2], [product_size; 2])?,
    )?;
    let dry = grinding
        .dry_throughput(
            feed_size,
            product_size,
            Power::new(1e6),
            Dimensionless::new(0.9),
        )?
        .value();
    let copper = ElementId::new("Cu")?;
    let mut classes = Vec::new();
    for (id, diameter, parts, grade) in [
        ("valuable-fine", 100e-6, 4.0, 0.3),
        ("gangue-fine", 100e-6, 6.0, 0.0),
        ("valuable-coarse", 200e-6, 6.0, 0.3),
        ("gangue-coarse", 200e-6, 14.0, 0.0),
    ] {
        let flow = dry * (parts / 30.0);
        classes.push(ParticleClass::new(
            id,
            Length::new(diameter),
            MassFlowRate::new(flow),
            BTreeMap::from([(copper.clone(), MassFlowRate::new(flow * grade))]),
        )?);
    }
    let feed = ParticleStream::new(classes, MassFlowRate::new(dry))?;
    let classified = classifier()?.classify(&feed, &QueryPoint::new())?;
    let conditions =
        QueryPoint::new().with_quantity("pH", QuantitySpec::dimensional(Dims::NONE), 10.0)?;
    let valuable = kinetics(true)?.evaluate(time, &conditions)?;
    let gangue = kinetics(false)?.evaluate(time, &conditions)?;
    let partitions = classified
        .reject
        .classes()
        .iter()
        .map(|class| {
            let response = if class.id().starts_with("valuable-") {
                valuable
            } else {
                gangue
            };
            Ok((
                class.id().to_owned(),
                Partition::from_fractions(
                    response.recovered,
                    Dimensionless::new(
                        response.recoverable_remaining.value()
                            + response.intrinsically_unrecoverable.value(),
                    ),
                )?,
            ))
        })
        .collect::<Result<Vec<_>, Box<dyn Error>>>()?;
    let recovered = classified.reject.split(
        &partitions,
        Partition::from_product(Dimensionless::new(0.1))?,
    )?;

    println!(
        "Authored single-pass scenario; PSD, assays and recovery rates are explicit illustrative inputs."
    );
    println!(
        "Bond energy={:.6} J/kg dry solids; net pinion power=1 MW; utilization=0.9.",
        grinding.specific_energy(feed_size, product_size)?.value()
    );
    println!(
        "Two equal steady mixed tanks; total residence={} s; pH=10; water recovery to concentrate=0.1.",
        time.value()
    );
    for (name, stream) in [
        ("Mill product", &feed),
        ("Classifier underflow", &classified.product),
        ("Classifier overflow / recovery feed", &classified.reject),
        ("Concentrate", &recovered.product),
        ("Tailings", &recovered.reject),
    ] {
        report(name, stream, &copper)?;
    }
    let outputs = [&classified.product, &recovered.product, &recovered.reject];
    let dry_out: f64 = outputs
        .iter()
        .map(|s| s.dry_solids().map(|q| q.value()))
        .collect::<Result<Vec<_>, _>>()?
        .iter()
        .sum();
    let cu_out: f64 = outputs
        .iter()
        .map(|s| s.element_flow(&copper).map(|q| q.value()))
        .collect::<Result<Vec<_>, _>>()?
        .iter()
        .sum();
    let water_out: f64 = outputs.iter().map(|s| s.water().value()).sum();
    println!(
        "Balance residuals [kg/s]: dry solids={:.3e}; Cu={:.3e}; water={:.3e}.",
        dry_out - feed.dry_solids()?.value(),
        cu_out - feed.element_flow(&copper)?.value(),
        water_out - feed.water().value()
    );
    println!(
        "Underflow remains a separate outlet; recycle convergence, liberation, entrainment and reagent response require their own models."
    );
    Ok(())
}
