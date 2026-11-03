//! Authored closed grinding/classification balances and a fresh-throughput Bond estimate.
//! Run `cargo run -p fs-material --example mineral_circuit -- --help` for usage.

use std::{collections::BTreeMap, error::Error};

use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance, QueryPoint,
};
use fs_material::comminution::{
    BondGrindingCard, BondSizeDomain, BondTestContext, BondWorkIndex,
    breakage::{
        BreakageColumn, BreakageTestContext, BreakageTransfer, SizeClass, SizeRedistributionCard,
    },
    circuit::{
        CircuitSettings, ClosedGrindingCircuit, MILL_DRY_FLOW_AXIS, power::BondCircuitPrototype,
    },
};
use fs_material::separation::{
    CLASSIFIER_LAW_ID, ClassifierCard, PARTICLE_DIAMETER_AXIS, ParticleClass, ParticleStream,
    Partition, SOLIDS_FRACTION_AXIS,
};
use fs_qty::{
    Dimensionless, Dims, Length, MassFlowRate, Power, QuantitySpec, Temperature,
    chemistry::ElementId,
};

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.len() == 1 && (args[0] == "--help" || args[0] == "-h") {
        println!(
            "Usage: mineral_circuit [coarse_retained_fraction [underflow_purge_fraction [net_pinion_kW]]]\n\
Defaults: coarse retention 0.5, underflow purge 0; both allowed in [0,1].\n\
The remaining coarse mass enters the fine class; unpurged underflow returns to the mill.\n\
Prints an authored steady recycle, conservative dry/Cu/water balances and discrete P80.\n\
The independent Bond estimate uses fresh throughput, F80 and combined overflow+purge P80.\n\
Optional positive net pinion kW also solves a scaled on-stream operating point and rechecks source limits.\n\
No fitted kinetics, predicted liberation or first-principles particle model is supplied."
        );
        return Ok(());
    }
    if args.len() > 3 {
        return Err(
            "expected optional coarse retention, underflow purge and net pinion kW; use --help"
                .into(),
        );
    }
    let retained: f64 = args.first().map_or(Ok(0.5), |value| value.parse())?;
    if !retained.is_finite() || !(0.0..=1.0).contains(&retained) {
        return Err("coarse retained fraction must be finite and in [0,1]".into());
    }
    let purge: f64 = args.get(1).map_or(Ok(0.0), |value| value.parse())?;
    let power_kw: Option<f64> = args.get(2).map(|value| value.parse()).transpose()?;
    if power_kw.is_some_and(|value| !value.is_finite() || value <= 0.0) {
        return Err("net pinion power must be finite and positive".into());
    }
    let underflow_recycle =
        Partition::from_fractions(Dimensionless::new(1.0 - purge), Dimensionless::new(purge))?;
    let description = format!(
        "Authored two-size scenario, not measured ore or a plant fit. Coarse retention={retained:.17e}; \
         Nonselective classifier-underflow purge={purge:.17e}, remainder returned to the mill. \
         Coarse 200 um: fresh 10 kg/s with Cu 2 kg/s; fine 100 um: fresh 2 kg/s with Cu 0.2 kg/s; water 12 kg/s. \
         Fine material retains its class; coarse remainder becomes fine, carrying the same assays. \
         Classifier d50 = 100 um, sharpness = 2, water bypass = 0.2. Constant per-pass coefficients are \
         explicitly assumed over actual mill dry flow 12..210 kg/s and solids mass fraction 0.5..0.95 \
         at 298.15 K. Separate Bond Wi = 15 kWh/dry metric tonne; all circuit rates are on-stream."
    );
    let artifact = hash_domain(
        "mineral-circuit-authored-scenario-v1",
        description.as_bytes(),
    );
    let provenance = Provenance {
        source: description.clone(),
        license: "project-authored example".into(),
        artifact: Some(artifact),
    };
    let quantity = QuantitySpec::dimensional;
    let common_validity = ValidityDomain::unconstrained()
        .with_quantity("temperature", quantity(Temperature::DIMS), 298.15, 298.15)
        .with_quantity(SOLIDS_FRACTION_AXIS, quantity(Dims::NONE), 0.5, 0.95);
    let mill = SizeRedistributionCard::new(
        [("coarse", 200e-6), ("fine", 100e-6)]
            .into_iter()
            .map(|(id, diameter)| SizeClass {
                id: id.into(),
                diameter: Length::new(diameter),
            })
            .collect(),
        [
            (
                "coarse",
                vec![("coarse", retained), ("fine", 1.0 - retained)],
            ),
            ("fine", vec![("fine", 1.0)]),
        ]
        .into_iter()
        .map(|(input, entries)| BreakageColumn {
            input: input.into(),
            transfers: entries
                .into_iter()
                .map(|(output, fraction)| BreakageTransfer {
                    output: output.into(),
                    mass_fraction: Dimensionless::new(fraction),
                })
                .collect(),
        })
        .collect(),
        BreakageTestContext {
            specimen: "authored two-class copper-bearing composite".into(),
            method: "authored constant congruent per-pass fractions; no fitted kinetics".into(),
            provenance: provenance.clone(),
        },
        common_validity.clone().with_quantity(
            MILL_DRY_FLOW_AXIS,
            quantity(MassFlowRate::DIMS),
            12.0,
            210.0,
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
        validity: common_validity.with_quantity(
            PARTICLE_DIAMETER_AXIS,
            quantity(Length::DIMS),
            100e-6,
            200e-6,
        ),
        sources: vec![artifact],
        provenance: provenance.clone(),
    })?;
    let copper = ElementId::new("Cu")?;
    let fresh = ParticleStream::new(
        [("coarse", 200e-6, 10.0, 2.0), ("fine", 100e-6, 2.0, 0.2)]
            .into_iter()
            .map(|(id, size, dry, cu)| {
                ParticleClass::new(
                    id,
                    Length::new(size),
                    MassFlowRate::new(dry),
                    BTreeMap::from([(copper.clone(), MassFlowRate::new(cu))]),
                )
            })
            .collect::<Result<Vec<_>, _>>()?,
        MassFlowRate::new(12.0),
    )?;
    let conditions =
        QueryPoint::new().with_quantity("temperature", quantity(Temperature::DIMS), 298.15)?;
    let circuit = ClosedGrindingCircuit::new(mill, classifier);
    let settings = CircuitSettings {
        absolute_balance_tolerance: MassFlowRate::new(1e-11),
        relative_balance_tolerance: Dimensionless::new(1e-12),
        maximum_circulating_load: Dimensionless::new(20.0),
    };
    let solved = circuit.solve_with_underflow_recycle(
        &fresh,
        &conditions,
        &conditions,
        settings,
        underflow_recycle,
    )?;
    let external = ParticleStream::mix(&[&solved.product, &solved.purge])?;

    println!("# {description}\n# Authored scenario source hash: {artifact}");
    println!(
        "# Unpurged underflow recycles. Overflow and purge are distinct external outlets. P80 is discrete, without within-bin interpolation."
    );
    println!("stream,dry_kg_s,Cu_kg_s,water_kg_s,P80_um");
    for (name, stream) in [
        ("fresh", &solved.fresh_feed),
        ("mill_feed", &solved.mill_feed),
        ("mill_discharge", &solved.mill_discharge),
        ("classifier_underflow", &solved.classifier_underflow),
        ("recycle", &solved.recycle),
        ("purge", &solved.purge),
        ("product", &solved.product),
        ("combined_external", &external),
    ] {
        let dry = stream.dry_solids()?.value();
        let passing = if dry == 0.0 {
            "undefined".to_owned()
        } else {
            format!(
                "{:.3}",
                1e6 * stream
                    .discrete_passing_size(Dimensionless::new(0.8))?
                    .value()
            )
        };
        println!(
            "{name},{dry:.9},{:.9},{:.9},{passing}",
            stream.element_flow(&copper)?.value(),
            stream.water().value(),
        );
    }
    println!(
        "Circulating load: {:.9} recycle/fresh dry mass ({:.3}%)",
        solved.circulating_load.value(),
        100.0 * solved.circulating_load.value()
    );
    println!(
        "Actual mill-feed solids mass fraction: {:.9}",
        solved.mill_feed.solids_mass_fraction()?.value()
    );
    println!(
        "Overflow + purge minus fresh residuals [kg/s]: dry={:.6e}, Cu={:.6e}, water={:.6e}",
        solved.closure.external_minus_fresh_solids.value(),
        solved.closure.external_minus_fresh_elements[&copper].value(),
        solved.closure.external_minus_fresh_water.value()
    );
    println!(
        "Maximum checked balance discrepancy: {:.6e} kg/s; relative={:.6e}",
        solved.closure.maximum_absolute_error.value(),
        solved.closure.maximum_relative_error.value()
    );

    let bond = BondGrindingCard::new(
        BondWorkIndex::from_kwh_per_metric_tonne(15.0)?,
        BondTestContext {
            specimen: "authored circuit work-index scenario".into(),
            method: "authored Wi and closing screen; no measured test or circuit calibration"
                .into(),
            closing_screen: Length::new(150e-6),
            provenance,
        },
        BondSizeDomain::new(
            [Length::new(100e-6), Length::new(200e-6)],
            [Length::new(100e-6), Length::new(200e-6)],
        )?,
    )?;
    let f80 = fresh.discrete_passing_size(Dimensionless::new(0.8))?;
    // Every dry outlet belongs to the circuit's mass/size boundary. Using only
    // the fine overflow would hide coarse material deliberately discharged.
    let p80 = external.discrete_passing_size(Dimensionless::new(0.8))?;
    if p80 < f80 {
        let energy = bond.specific_energy(f80, p80)?.value();
        // Circuit-scale Bond energy uses fresh throughput, not circulating mill feed.
        let required_pinion_power = energy * fresh.dry_solids()?.value();
        println!(
            "Independent Bond estimate to combined overflow+purge P80: {:.6} J/kg dry ({:.6} kWh/dry metric tonne)",
            energy,
            energy / 3600.0
        );
        println!(
            "Required net pinion power: {:.6} kW at on-stream fresh {:.6} kg/s",
            required_pinion_power / 1000.0,
            fresh.dry_solids()?.value()
        );
    } else {
        println!(
            "Bond estimate unavailable: combined external P80 equals fresh F80; this law requires strict size reduction."
        );
    }
    if let Some(power_kw) = power_kw {
        let prototype = BondCircuitPrototype::new(
            &circuit,
            &bond,
            &fresh,
            &conditions,
            &conditions,
            settings,
            underflow_recycle,
        )?;
        let point = prototype.at_net_pinion_power(Power::new(power_kw * 1000.0))?;
        println!(
            "On-stream Bond operating point: {:.6} kW; F80={:.3} um, combined external P80={:.3} um, energy={:.6} J/kg dry",
            point.net_pinion_power.value() / 1000.0,
            point.fresh_f80.value() * 1e6,
            point.external_p80.value() * 1e6,
            point.specific_energy.value(),
        );
        println!("scaled_stream,dry_kg_s,Cu_kg_s,water_kg_s");
        for (name, stream) in [
            ("fresh", &point.solution.fresh_feed),
            ("mill_feed", &point.solution.mill_feed),
            ("recycle", &point.solution.recycle),
            ("overflow", &point.solution.product),
            ("purge", &point.solution.purge),
        ] {
            println!(
                "{name},{:.9},{:.9},{:.9}",
                stream.dry_solids()?.value(),
                stream.element_flow(&copper)?.value(),
                stream.water().value(),
            );
        }
        println!(
            "Rechecked actual mill flow {:.6} kg/s, solids fraction {:.9}, circulating load {:.9}; maximum balance discrepancy {:.6e} kg/s",
            point.solution.mill_feed.dry_solids()?.value(),
            point.solution.mill_feed.solids_mass_fraction()?.value(),
            point.solution.circulating_load.value(),
            point.solution.closure.maximum_absolute_error.value(),
        );
        println!(
            "Scaled rates are on-stream; reporting-period utilization is separate from source operating limits."
        );
    }
    println!(
        "Per-pass fractions and Wi are independent authored inputs; no kinetics, liberation or operating-efficiency prediction."
    );
    Ok(())
}
