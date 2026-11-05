//! cargo run -p fs-material --example transient_grinding_circuit
//!
//! Synthetic feed-grade changes through a retained mill/classifier inventory.
//! All rates, partitions and support ranges below are authored, not measured.
//! Instantaneous recycle has no sump/pipe hold-up; the existing local grinding
//! and classifier owners are used at every first-order numerical time step.

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
            ContinuousPopulationBalanceMill, MEAN_RESIDENCE_TIME_AXIS, MILL_DRY_FLOW_AXIS,
            SelectionBreakageClass,
            circuit::{ContinuousGrindingCircuit, transient::TransientCircuitSettings},
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

fn feed(copper: &ElementId, copper_rate: f64) -> Result<ParticleStream, Box<dyn Error>> {
    Ok(ParticleStream::new(
        [
            ("coarse", 200e-6, 10.0, copper_rate),
            ("fine", 100e-6, 0.0, 0.0),
        ]
        .into_iter()
        .map(|(id, diameter, dry, cu)| {
            ParticleClass::new(
                id,
                Length::new(diameter),
                MassFlowRate::new(dry),
                BTreeMap::from([(copper.clone(), MassFlowRate::new(cu))]),
            )
        })
        .collect::<Result<_, _>>()?,
        MassFlowRate::new(5.0),
    )?)
}

fn main() -> Result<(), Box<dyn Error>> {
    let artifact = hash_bytes(b"authored-transient-grinding-example-v1");
    let source = Provenance {
        source: "authored two-class transient example; no ore or plant calibration".into(),
        license: "project example".into(),
        artifact: Some(artifact),
    };
    let mill = ContinuousPopulationBalanceMill::new(
        [("coarse", 200e-6), ("fine", 100e-6)]
            .into_iter()
            .map(|(id, size)| SizeClass {
                id: id.into(),
                diameter: Length::new(size),
            })
            .collect(),
        vec![
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
        ],
        BreakageTestContext {
            specimen: "synthetic two-class feed".into(),
            method: "authored continuous selection/daughter parameters".into(),
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
        .map(|(id, value, dims)| (id.into(), LawParameter { value, dims }))
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
                0.0,
                1.0,
            ),
        sources: vec![artifact],
        provenance: source,
    })?;
    let circuit = ContinuousGrindingCircuit::new(mill, classifier);
    let copper = ElementId::new("Cu")?;
    let baseline = feed(&copper, 1.0)?;
    let residence = Time::new(10.0);
    let returned = Partition::from_product(Dimensionless::new(0.75))?;
    let conditions = QueryPoint::new();
    let steady = circuit.at_residence_time(
        &baseline,
        residence,
        &conditions,
        &conditions,
        CircuitSettings {
            absolute_balance_tolerance: MassFlowRate::new(1e-9),
            relative_balance_tolerance: Dimensionless::new(1e-11),
            maximum_circulating_load: Dimensionless::new(10.0),
        },
        returned,
    )?;
    let mut inventory = MillInventory::new(steady.mill.class_holdup, steady.mill.water_holdup)?;
    let settings = TransientCircuitSettings {
        time: TransientMillSettings {
            maximum_step: Time::new(0.25),
            maximum_steps: 1000,
            absolute_balance_tolerance: Mass::new(1e-8),
            relative_balance_tolerance: Dimensionless::new(1e-11),
        },
        maximum_circulating_load: Dimensionless::new(10.0),
    };
    println!(
        "SYNTHETIC closed grinding transient: fixed kinetics, instantaneous classifier recycle."
    );
    println!("Fresh dry flow 10 kg/s; water 5 kg/s; mill residence 10 s; underflow returned 75%.");
    println!(
        "Initial retained inventory: {:.6} kg dry, {:.6} kg H2O, {:.6} kg Cu.",
        inventory.dry_solids()?.value(),
        inventory.water().value(),
        inventory.element_mass(&copper)?.value()
    );
    for (label, cu, seconds) in [
        ("first second after Cu doubles", 2.0, 1.0),
        ("next 49 seconds at higher Cu", 2.0, 49.0),
        ("ten seconds after Cu returns", 1.0, 10.0),
    ] {
        let fresh = feed(&copper, cu)?;
        let step = circuit.advance_inventory(
            &inventory,
            &fresh,
            residence,
            Time::new(seconds),
            &conditions,
            &conditions,
            returned,
            settings,
            || false,
        )?;
        let incoming = step.integrated_fresh.element_mass(&copper)?.value();
        let external = step.integrated_product.element_mass(&copper)?.value()
            + step.integrated_purge.element_mass(&copper)?.value();
        let accumulation = step.inventory.element_mass(&copper)?.value()
            - inventory.element_mass(&copper)?.value();
        println!(
            "\n{label}: {seconds:.2} s, {} accepted steps",
            step.substeps
        );
        println!(
            "  Cu in={incoming:.9} kg; overflow+purge={external:.9} kg; retained change={accumulation:.9} kg"
        );
        println!(
            "  Same-period output/input={:.6}; boundary Cu residual={:.3e} kg",
            external / incoming,
            step.boundary_balance.elements[&copper].value()
        );
        println!(
            "  Final overflow Cu={:.9} kg/s; dry circulating load={:.6} return/fresh",
            step.final_point.product.element_flow(&copper)?.value(),
            step.final_point.circulating_load.value()
        );
        inventory = step.inventory;
    }
    println!(
        "\nChanging same-period output/input is inventory transport, not inferred metal loss or recovery."
    );
    println!(
        "This is a bounded first-order trajectory, not an error enclosure or a plant-calibration claim."
    );
    Ok(())
}
