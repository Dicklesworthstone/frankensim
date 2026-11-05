//! cargo run -p fs-material --example grinding_sump
//!
//! Authored two-class wet grinding with a finite classifier-feed sump. No ore,
//! source rate, tank geometry or industrial operating calibration is inferred.
//! All three intervals succeed before any successful result is printed.

use fs_blake3::hash_bytes;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance, QueryPoint,
};
use fs_material::{
    comminution::{
        breakage::{BreakageTestContext, BreakageTransfer, SizeClass},
        circuit::CircuitSettings,
        population_balance::{
            ClassHoldup, ContinuousPopulationBalanceMill, MEAN_RESIDENCE_TIME_AXIS,
            MILL_DRY_FLOW_AXIS, SelectionBreakageClass,
            circuit::{
                ContinuousGrindingCircuit,
                transient::{
                    TransientCircuitSettings,
                    sump::{ClassifierFeedSump, SumpCircuitInventory},
                },
            },
            transient::{MillInventory, TransientMillSettings},
        },
    },
    separation::{
        CLASSIFIER_LAW_ID, ClassifierCard, PARTICLE_DIAMETER_AXIS, ParticleClass, ParticleStream,
        Partition, SOLIDS_FRACTION_AXIS,
    },
};
use fs_qty::{
    Dimensionless, Dims, Frequency, Length, Mass, MassFlowRate, QuantitySpec, Time,
    chemistry::ElementId,
};
use std::{collections::BTreeMap, error::Error};

fn feed(copper: &ElementId, cu: f64) -> Result<ParticleStream, Box<dyn Error>> {
    Ok(ParticleStream::new(
        [("coarse", 200e-6, 10.0, cu), ("fine", 100e-6, 0.0, 0.0)]
            .into_iter()
            .map(|(id, size, dry, copper_flow)| {
                ParticleClass::new(
                    id,
                    Length::new(size),
                    MassFlowRate::new(dry),
                    BTreeMap::from([(copper.clone(), MassFlowRate::new(copper_flow))]),
                )
            })
            .collect::<Result<_, _>>()?,
        MassFlowRate::new(5.0),
    )?)
}

fn main() -> Result<(), Box<dyn Error>> {
    let unsupported = match std::env::args().skip(1).collect::<Vec<_>>().as_slice() {
        [] => false,
        [flag] if flag == "--unsupported-dilution" => true,
        _ => return Err("usage: grinding_sump [--unsupported-dilution]".into()),
    };
    let source = Provenance {
        source: "authored mill/classifier/sump illustration; not plant testwork".into(),
        license: "project example".into(),
        artifact: Some(hash_bytes(b"authored-finite-grinding-sump-v1")),
    };
    let classes = [("coarse", 200e-6), ("fine", 100e-6)]
        .into_iter()
        .map(|(id, d)| SizeClass {
            id: id.into(),
            diameter: Length::new(d),
        })
        .collect();
    let kinetics = vec![
        SelectionBreakageClass {
            input: "coarse".into(),
            selection_rate: Frequency::new(0.1),
            daughters: vec![BreakageTransfer {
                output: "fine".into(),
                mass_fraction: Dimensionless::new(1.0),
            }],
        },
        SelectionBreakageClass {
            input: "fine".into(),
            selection_rate: Frequency::new(0.0),
            daughters: vec![],
        },
    ];
    let mill = ContinuousPopulationBalanceMill::new(
        classes,
        kinetics,
        BreakageTestContext {
            specimen: "synthetic two-size ore".into(),
            method: "authored selection/daughter parameters".into(),
            provenance: source.clone(),
        },
        ValidityDomain::unconstrained()
            .with_quantity(
                MEAN_RESIDENCE_TIME_AXIS,
                QuantitySpec::dimensional(Time::DIMS),
                1.0,
                100.0,
            )
            .with_quantity(
                MILL_DRY_FLOW_AXIS,
                QuantitySpec::dimensional(MassFlowRate::DIMS),
                1.0,
                1000.0,
            )
            .with_quantity(
                SOLIDS_FRACTION_AXIS,
                QuantitySpec::dimensional(Dims::NONE),
                0.0,
                1.0,
            ),
    )?;
    let classifier = ClassifierCard::from_model_card(&ConstitutiveModelCard {
        law: LawId(CLASSIFIER_LAW_ID.into()),
        law_version: 1,
        state_schema_version: 1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        parameters: [
            ("corrected-cut-size", 100e-6, Length::DIMS),
            ("sharpness", 1.0, Dims::NONE),
            ("water-recovery-to-underflow", 0.2, Dims::NONE),
        ]
        .into_iter()
        .map(|(name, value, dims)| (name.into(), LawParameter { value, dims }))
        .collect(),
        validity: ValidityDomain::unconstrained()
            .with_quantity(
                PARTICLE_DIAMETER_AXIS,
                QuantitySpec::dimensional(Length::DIMS),
                100e-6,
                200e-6,
            )
            .with_quantity(
                SOLIDS_FRACTION_AXIS,
                QuantitySpec::dimensional(Dims::NONE),
                0.5,
                0.9,
            ),
        sources: vec![source.artifact.unwrap()],
        provenance: source.clone(),
    })?;
    let model = ContinuousGrindingCircuit::new(mill, classifier);
    let sump = ClassifierFeedSump::new(Time::new(4.0), source)?;
    let cu = ElementId::new("Cu")?;
    let residence = Time::new(10.0);
    let returned = Partition::from_product(Dimensionless::new(0.75))?;
    let query = QueryPoint::new();
    let steady = model.at_residence_time(
        &feed(&cu, 1.0)?,
        residence,
        &query,
        &query,
        CircuitSettings {
            absolute_balance_tolerance: MassFlowRate::new(1e-9),
            relative_balance_tolerance: Dimensionless::new(1e-11),
            maximum_circulating_load: Dimensionless::new(10.0),
        },
        returned,
    )?;
    // The steady sump stock is the original stream times its explicit drain time.
    let sump_stock = MillInventory::new(
        steady
            .mill
            .outlet
            .classes()
            .iter()
            .map(|c| ClassHoldup {
                id: c.id().into(),
                diameter: c.diameter(),
                dry_solids: Mass::new(c.dry_solids().value() * sump.residence().value()),
                elements: c
                    .elements()
                    .iter()
                    .map(|(e, m)| (e.clone(), Mass::new(m.value() * sump.residence().value())))
                    .collect(),
            })
            .collect(),
        Mass::new(steady.mill.outlet.water().value() * sump.residence().value()),
    )?;
    let mut inventory = SumpCircuitInventory {
        mill: MillInventory::new(steady.mill.class_holdup, steady.mill.water_holdup)?,
        sump: sump_stock,
    };
    let settings = TransientCircuitSettings {
        time: TransientMillSettings {
            maximum_step: Time::new(0.25),
            maximum_steps: 1000,
            absolute_balance_tolerance: Mass::new(1e-8),
            relative_balance_tolerance: Dimensionless::new(1e-11),
        },
        maximum_circulating_load: Dimensionless::new(10.0),
    };
    let mut results = Vec::new();
    for (label, grade, makeup, seconds) in [
        (
            "first second after Cu doubles and sump dilution starts",
            2.0,
            3.0,
            1.0,
        ),
        ("next 29 seconds under those conditions", 2.0, 3.0, 29.0),
        (
            "15 seconds after original Cu feed and no sump dilution resume",
            1.0,
            0.0,
            15.0,
        ),
    ] {
        let makeup = if unsupported { 1e7 } else { makeup };
        let step = model.advance_with_sump(
            &sump,
            &inventory,
            &feed(&cu, grade)?,
            MassFlowRate::new(makeup),
            residence,
            Time::new(seconds),
            &query,
            &query,
            returned,
            settings,
            || false,
        )?;
        inventory = step.inventory.clone();
        results.push((label, step));
    }
    println!(
        "SYNTHETIC wet grinding: mill residence 10 s, finite sump residence 4 s, underflow return 75%."
    );
    println!(
        "Source rates/partitions and sump drain law are explicit assumptions, not plant calibration."
    );
    for (label, step) in results {
        let in_cu = step.integrated_fresh.element_mass(&cu)?.value();
        let out_cu = step.integrated_product.element_mass(&cu)?.value()
            + step.integrated_purge.element_mass(&cu)?.value();
        let retained = step.inventory.mill.element_mass(&cu)?.value()
            + step.inventory.sump.element_mass(&cu)?.value();
        let initial = step.initial().mill.element_mass(&cu)?.value()
            + step.initial().sump.element_mass(&cu)?.value();
        println!(
            "\n{label}: {:.2} s, {} steps",
            step.elapsed.value(),
            step.substeps
        );
        println!(
            "  Cu in={in_cu:.9} kg; external={out_cu:.9} kg; stock change={:.9} kg",
            retained - initial
        );
        println!(
            "  Sump dilution={:.6} kg; mill water stock={:.6} kg; sump water stock={:.6} kg",
            step.integrated_makeup_water.value(),
            step.inventory.mill.water().value(),
            step.inventory.sump.water().value()
        );
        println!(
            "  Final classifier solids fraction={:.6}; mill-discharge Cu={:.9} kg/s; classifier-feed Cu={:.9} kg/s",
            step.final_point
                .classifier_feed
                .solids_mass_fraction()?
                .value(),
            step.final_point.mill_discharge.element_flow(&cu)?.value(),
            step.final_point.classifier_feed.element_flow(&cu)?.value()
        );
        println!(
            "  Cu residuals: mill={:.3e} kg, sump={:.3e} kg, whole={:.3e} kg",
            step.mill_balance.elements[&cu].value(),
            step.sump_balance.elements[&cu].value(),
            step.boundary_balance.elements[&cu].value()
        );
    }
    println!(
        "\nFirst-order endpoint-checked trajectory; no slip, pump/level law or continuous-time enclosure."
    );
    Ok(())
}
