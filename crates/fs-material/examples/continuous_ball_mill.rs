//! Throughput versus product PSD for a source-declared continuous grinding law.
//!
//! cargo run -p fs-material --example continuous_ball_mill
//!
//! Every numeric rate, daughter fraction, operating range and feed in this
//! example is authored synthetic data. Nothing is a measured ore constant or
//! a calibrated plant prediction. Hold-up is prescribed, not inferred from
//! geometry or power. Water shares the solids' ideal mixed residence time.

use fs_evidence::ValidityDomain;
use fs_matdb::{Provenance, QueryPoint};
use fs_material::{
    comminution::{
        breakage::{BreakageTestContext, BreakageTransfer, SizeClass},
        population_balance::{
            ContinuousPopulationBalanceMill, MEAN_RESIDENCE_TIME_AXIS, MILL_DRY_FLOW_AXIS,
            SelectionBreakageClass,
        },
    },
    separation::{ParticleClass, ParticleStream, SOLIDS_FRACTION_AXIS},
};
use fs_qty::{
    Dimensionless, Frequency, Length, Mass, MassFlowRate, QuantitySpec, Time, chemistry::ElementId,
};
use std::{collections::BTreeMap, error::Error};

fn main() -> Result<(), Box<dyn Error>> {
    let layout: Vec<_> = [("coarse", 1000e-6), ("medium", 250e-6), ("fine", 75e-6)]
        .into_iter()
        .map(|(id, size)| SizeClass {
            id: id.into(),
            diameter: Length::new(size),
        })
        .collect();
    let daughter = |output: &str, fraction| BreakageTransfer {
        output: output.into(),
        mass_fraction: Dimensionless::new(fraction),
    };
    let model = ContinuousPopulationBalanceMill::new(
        layout.clone(),
        vec![
            SelectionBreakageClass { input: "coarse".into(), selection_rate: Frequency::new(0.02),
                daughters: vec![daughter("medium", 0.7), daughter("fine", 0.3)] },
            SelectionBreakageClass { input: "medium".into(), selection_rate: Frequency::new(0.05),
                daughters: vec![daughter("fine", 1.0)] },
            SelectionBreakageClass { input: "fine".into(), selection_rate: Frequency::new(0.0), daughters: vec![] },
        ],
        BreakageTestContext {
            specimen: "synthetic ore for continuous ball-mill illustration".into(),
            method: "authored local selection rates and daughter mass fractions; no measured calibration".into(),
            provenance: Provenance { source: "continuous_ball_mill authored demonstration".into(),
                license: "project example".into(), artifact: None },
        },
        ValidityDomain::unconstrained()
            .with_quantity(MEAN_RESIDENCE_TIME_AXIS, QuantitySpec::dimensional(Time::DIMS), 10.0, 120.0)
            .with_quantity(MILL_DRY_FLOW_AXIS, QuantitySpec::dimensional(MassFlowRate::DIMS), 50.0, 300.0)
            .with_quantity(SOLIDS_FRACTION_AXIS, QuantitySpec::dimensional(Dimensionless::DIMS), 0.5, 0.9),
    )?;
    let copper = ElementId::new("Cu")?;
    println!("SYNTHETIC continuous ball-mill example; no measured ore or plant calibration.");
    println!(
        "Fixed dry hold-up: 6000 kg; source rates and daughters fixed within their declared support."
    );
    println!("dry feed kg/s | residence s | coarse % | medium % | fine % | Cu out kg/s");
    for throughput in [100.0, 200.0] {
        let classes = layout
            .iter()
            .zip([0.8, 0.15, 0.05])
            .map(|(size, fraction)| {
                let dry = throughput * fraction;
                ParticleClass::new(
                    size.id.clone(),
                    size.diameter,
                    MassFlowRate::new(dry),
                    BTreeMap::from([(copper.clone(), MassFlowRate::new(dry * 0.01))]),
                )
            })
            .collect::<Result<Vec<_>, _>>()?;
        let feed = ParticleStream::new(classes, MassFlowRate::new(throughput * 0.5))?;
        let solution = model.at_dry_holdup(&feed, Mass::new(6000.0), &QueryPoint::new())?;
        let percent = |name: &str| {
            solution
                .outlet
                .classes()
                .iter()
                .find(|c| c.id() == name)
                .expect("source class is retained")
                .dry_solids()
                .value()
                / throughput
                * 100.0
        };
        println!(
            "{throughput:13.1} | {:11.1} | {:8.2} | {:8.2} | {:6.2} | {:11.3}",
            solution.mean_residence_time.value(),
            percent("coarse"),
            percent("medium"),
            percent("fine"),
            solution.outlet.element_flow(&copper)?.value()
        );
        println!(
            "  balanced outlet: {:.1} kg/s dry + {:.1} kg/s water; hold-up: {:.1} kg dry + {:.1} kg water",
            solution.outlet.dry_solids()?.value(),
            solution.outlet.water().value(),
            solution.dry_solids_holdup.value(),
            solution.water_holdup.value()
        );
    }
    println!(
        "At fixed hold-up, the faster feed has less grinding residence and a coarser product."
    );
    Ok(())
}
