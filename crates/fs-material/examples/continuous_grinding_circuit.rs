//! cargo run -p fs-material --example continuous_grinding_circuit
//!
//! Authored continuous mill/classifier example, not measured ore or plant data.
//! Selection frequencies, daughters, classifier coefficients and applicability
//! ranges are explicitly synthetic. Prescribed residence reports its resulting
//! inventory; prescribed dry hold-up solves residence against circulating feed
//! within an explicit valid bracket. Neither is inferred from a Bond work index.
//! A separate authored Bond card estimates energy at each actual operating PSD.

use fs_blake3::hash_bytes;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance, QueryPoint,
};
use fs_material::{
    comminution::{
        BondGrindingCard, BondSizeDomain, BondTestContext, BondWorkIndex,
        breakage::{BreakageTestContext, BreakageTransfer, SizeClass},
        circuit::CircuitSettings,
        population_balance::{
            ContinuousPopulationBalanceMill, MEAN_RESIDENCE_TIME_AXIS, MILL_DRY_FLOW_AXIS,
            SelectionBreakageClass,
            circuit::{
                ContinuousGrindingCircuit, holdup::HoldupSolveSettings,
                power::ContinuousBondCircuit,
            },
        },
    },
    separation::{
        CLASSIFIER_LAW_ID, ClassifierCard, PARTICLE_DIAMETER_AXIS, ParticleClass, ParticleStream,
        Partition, SOLIDS_FRACTION_AXIS,
    },
};
use fs_qty::{
    Dimensionless, Dims, Frequency, Length, Mass, MassFlowRate, Pressure, QuantitySpec, Time,
    chemistry::ElementId,
};
use std::{collections::BTreeMap, error::Error};

fn main() -> Result<(), Box<dyn Error>> {
    let source = Provenance {
        source: "authored continuous grinding/recycle demonstration; no measured calibration"
            .into(),
        license: "project example".into(),
        artifact: Some(hash_bytes(b"continuous-grinding-circuit-example-v1")),
    };
    let mill = ContinuousPopulationBalanceMill::new(
        vec![
            SizeClass {
                id: "coarse".into(),
                diameter: Length::new(200e-6),
            },
            SizeClass {
                id: "fine".into(),
                diameter: Length::new(100e-6),
            },
        ],
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
            specimen: "synthetic two-class assayed feed".into(),
            method: "authored continuous selection frequency and daughter fractions".into(),
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
            ("sharpness", 2.0, Dims::NONE),
            ("water-recovery-to-underflow", 0.2, Dims::NONE),
        ]
        .into_iter()
        .map(|(name, value, dims)| (name.into(), LawParameter { value, dims }))
        .collect(),
        validity: ValidityDomain::unconstrained()
            .with_quantity(
                PARTICLE_DIAMETER_AXIS,
                QuantitySpec::dimensional(Length::DIMS),
                50e-6,
                1e-3,
            )
            .with_quantity(
                SOLIDS_FRACTION_AXIS,
                QuantitySpec::dimensional(Dims::NONE),
                0.0,
                1.0,
            )
            .with_quantity(
                "pressure",
                QuantitySpec::dimensional(Pressure::DIMS),
                100_000.0,
                100_000.0,
            ),
        sources: vec![source.artifact.unwrap()],
        provenance: source,
    })?;
    let copper = ElementId::new("Cu")?;
    let fresh = ParticleStream::new(
        vec![
            ParticleClass::new(
                "coarse",
                Length::new(200e-6),
                MassFlowRate::new(10.0),
                BTreeMap::from([(copper.clone(), MassFlowRate::new(2.0))]),
            )?,
            ParticleClass::new(
                "fine",
                Length::new(100e-6),
                MassFlowRate::new(0.0),
                BTreeMap::from([(copper.clone(), MassFlowRate::new(0.0))]),
            )?,
        ],
        MassFlowRate::new(12.0),
    )?;
    let circuit = ContinuousGrindingCircuit::new(mill, classifier);
    let classifier_conditions = QueryPoint::new().with_quantity(
        "pressure",
        QuantitySpec::dimensional(Pressure::DIMS),
        100_000.0,
    )?;
    println!(
        "SYNTHETIC continuous grinding/classifier circuit; no measured ore or plant calibration."
    );
    println!("Fresh feed: 10 kg/s coarse dry solids, 2 kg/s assayed Cu, 12 kg/s water.");
    println!("The same continuous selection/daughter law is used at every prescribed residence.");
    let balance_settings = CircuitSettings {
        absolute_balance_tolerance: MassFlowRate::new(1e-10),
        relative_balance_tolerance: Dimensionless::new(1e-12),
        maximum_circulating_load: Dimensionless::new(100.0),
    };
    for (seconds, returned) in [(5.0, 1.0), (10.0, 1.0), (20.0, 1.0), (10.0, 0.75)] {
        let solution = circuit.at_residence_time(
            &fresh,
            Time::new(seconds),
            &QueryPoint::new(),
            &classifier_conditions,
            balance_settings,
            Partition::from_product(Dimensionless::new(returned))?,
        )?;
        let flows = &solution.circuit;
        let fine = flows
            .product
            .classes()
            .iter()
            .find(|class| class.id() == "fine")
            .expect("declared fine class");
        println!("\nResidence {seconds:.1} s; underflow returned {returned:.2}:");
        println!(
            "  Mill dry feed {:.9} kg/s; circulating load {:.9} recycle/fresh",
            flows.mill_feed.dry_solids()?.value(),
            flows.circulating_load.value()
        );
        println!(
            "  Overflow {:.9} kg/s dry, including {:.9} kg/s fine",
            flows.product.dry_solids()?.value(),
            fine.dry_solids().value()
        );
        println!(
            "  Purge {:.9} kg/s dry; overflow Cu {:.9} kg/s; purge Cu {:.9} kg/s",
            flows.purge.dry_solids()?.value(),
            flows.product.element_flow(&copper)?.value(),
            flows.purge.element_flow(&copper)?.value()
        );
        println!(
            "  Water: mill {:.9}, recycle {:.9}, overflow {:.9}, purge {:.9} kg/s",
            flows.mill_feed.water().value(),
            flows.recycle.water().value(),
            flows.product.water().value(),
            flows.purge.water().value()
        );
        println!(
            "  Predicted hold-up: {:.9} kg dry, {:.9} kg water",
            solution.mill.dry_solids_holdup.value(),
            solution.mill.water_holdup.value()
        );
        println!(
            "  Maximum checked account residual {:.3e} kg/s",
            flows.closure.maximum_absolute_error.value()
        );
    }
    println!(
        "\nFixed dry hold-up 333.333333333 kg; full underflow recycle; residence bracket 1..20 s:"
    );
    let powered = ContinuousBondCircuit::new(
        circuit.clone(),
        BondGrindingCard::new(
            BondWorkIndex::from_kwh_per_metric_tonne(15.0)?,
            BondTestContext {
                specimen: "synthetic two-class assayed feed".into(),
                method: "independently authored Bond work-index example".into(),
                closing_screen: Length::new(125e-6),
                provenance: Provenance {
                    source: "authored 15 kWh/metric-tonne Bond input; not measured testwork".into(),
                    license: "project example".into(),
                    artifact: Some(hash_bytes(b"continuous-grinding-bond-example-v1")),
                },
            },
            BondSizeDomain::new(
                [Length::new(200e-6), Length::new(1e-3)],
                [Length::new(70e-6), Length::new(150e-6)],
            )?,
        )?,
    );
    println!(
        "Separate authored Bond input: Wi 15 kWh/metric-tonne; actual external P80 is recomputed."
    );
    let doubled = ParticleStream::mix(&[&fresh, &fresh])?;
    for feed in [&fresh, &doubled] {
        let point = powered.at_dry_holdup(
            feed,
            Mass::new(333.333333333),
            &QueryPoint::new(),
            &classifier_conditions,
            balance_settings,
            Partition::from_product(Dimensionless::new(1.0))?,
            HoldupSolveSettings {
                minimum_residence: Time::new(1.0),
                maximum_residence: Time::new(20.0),
                absolute_mass_tolerance: Mass::new(1e-9),
                relative_mass_tolerance: Dimensionless::new(1e-12),
                maximum_iterations: 80,
            },
        )?;
        let fixed = point
            .holdup_solution()
            .expect("requested fixed dry hold-up");
        let flows = &fixed.solution.circuit;
        let fine = flows
            .product
            .classes()
            .iter()
            .find(|class| class.id() == "fine")
            .expect("declared fine class");
        println!(
            "  Fresh {:.1} kg/s dry -> residence {:.9} s, mill dry feed {:.9} kg/s, overflow {:.6}% fine",
            feed.dry_solids()?.value(),
            fixed.solution.mill.mean_residence_time.value(),
            flows.mill_feed.dry_solids()?.value(),
            100.0 * fine.dry_solids().value() / flows.product.dry_solids()?.value(),
        );
        println!(
            "    Original dry hold-up {:.9} kg; water hold-up {:.9} kg; dry residual {:.3e} kg after {} midpoint trials",
            fixed.solution.mill.dry_solids_holdup.value(),
            fixed.solution.mill.water_holdup.value(),
            fixed.mass_residual.value(),
            fixed.iterations,
        );
        println!(
            "    F80 {:.1} um; combined overflow/purge P80 {:.1} um; Bond {:.6} kWh/t; estimated net pinion {:.6} kW",
            point.fresh_f80.value() * 1e6,
            point.external_p80.value() * 1e6,
            point.specific_energy.value() / 3600.0,
            point.net_pinion_power.value() / 1000.0,
        );
    }
    println!(
        "Purge is a nonselective outlet. Bond energy does not calibrate selection rates or infer equipment power draw, liberation or reactive mineralogy."
    );
    Ok(())
}
