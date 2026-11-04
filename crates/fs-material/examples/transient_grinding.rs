//! Startup, feed change and washout of one perfectly mixed grinding volume.
//!
//! Every rate and operating window here is an authored illustration, not
//! specimen testwork or a calibrated plant. No Bond work index is converted
//! into a breakage frequency or transient power draw.

use fs_evidence::ValidityDomain;
use fs_matdb::{Provenance, QueryPoint};
use fs_material::{
    comminution::{
        breakage::{BreakageTestContext, BreakageTransfer, SizeClass},
        population_balance::{
            ContinuousPopulationBalanceMill, MEAN_RESIDENCE_TIME_AXIS, MILL_DRY_FLOW_AXIS,
            SelectionBreakageClass,
            transient::{MillInventory, TransientMillSettings},
        },
    },
    separation::{ParticleClass, ParticleStream},
};
use fs_qty::{
    Dimensionless, Frequency, Length, Mass, MassFlowRate, QuantitySpec, Time, chemistry::ElementId,
};
use std::collections::BTreeMap;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let copper = ElementId::new("Cu")?;
    let layout = [("coarse", 300e-6, 8.0, 0.02), ("fine", 75e-6, 2.0, 0.01)];
    let model = ContinuousPopulationBalanceMill::new(
        layout
            .iter()
            .map(|&(id, diameter, _, _)| SizeClass {
                id: id.into(),
                diameter: Length::new(diameter),
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
            specimen: "synthetic assayed two-size feed; not measured ore".into(),
            method: "authored constant selection and daughter coefficients".into(),
            provenance: Provenance {
                source: "analytical perfectly mixed inventory illustration".into(),
                license: "project fixture".into(),
                artifact: None,
            },
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
                0.0,
                100.0,
            ),
    )?;
    let feed = ParticleStream::new(
        layout
            .iter()
            .map(|&(id, diameter, flow, grade)| {
                ParticleClass::new(
                    id,
                    Length::new(diameter),
                    MassFlowRate::new(flow),
                    BTreeMap::from([(copper.clone(), MassFlowRate::new(flow * grade))]),
                )
            })
            .collect::<Result<_, _>>()?,
        MassFlowRate::new(5.0),
    )?;
    let settings = TransientMillSettings {
        maximum_step: Time::new(0.05),
        maximum_steps: 1000,
        absolute_balance_tolerance: Mass::new(1e-9),
        relative_balance_tolerance: Dimensionless::new(1e-12),
    };
    println!("Synthetic local grinding response; no measured ore/plant calibration.");
    println!("Common residence 10 s; coarse selection 0.1 /s; all daughters enter fines.");
    println!(
        "The declared source permits zero-flow washout and does not require an undefined empty-feed solids fraction."
    );
    println!(
        "phase,time_s,feed_dry_kg_s,held_dry_kg,held_coarse_kg,held_fine_kg,held_water_kg,held_cu_kg,discharged_dry_kg,dry_residual_kg"
    );
    let mut inventory = MillInventory::empty_like(&feed)?;
    let mut fed = 0.0;
    let mut discharged = 0.0;
    for (index, (phase, factor)) in [("startup", 1.0), ("double_feed", 2.0), ("washout", 0.0)]
        .into_iter()
        .enumerate()
    {
        let incoming = feed.scaled(Dimensionless::new(factor))?;
        let step = model.advance_inventory(
            &inventory,
            &incoming,
            Time::new(10.0),
            Time::new(30.0),
            &QueryPoint::new(),
            settings,
            || false,
        )?;
        let class_mass = |id: &str| {
            step.inventory
                .classes()
                .iter()
                .find(|c| c.id == id)
                .expect("original two-class layout")
                .dry_solids
                .value()
        };
        fed += step.integrated_feed.dry_solids()?.value();
        discharged += step.integrated_discharge.dry_solids()?.value();
        println!(
            "{phase},{},{:.6},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.3e}",
            30 * (index + 1),
            incoming.dry_solids()?.value(),
            step.inventory.dry_solids()?.value(),
            class_mass("coarse"),
            class_mass("fine"),
            step.inventory.water().value(),
            step.inventory.element_mass(&copper)?.value(),
            step.integrated_discharge.dry_solids()?.value(),
            step.balance.dry_solids.value(),
        );
        inventory = step.inventory;
    }
    println!(
        "All intervals: fed={fed:.9} kg, discharged={discharged:.9} kg, retained={:.9} kg, residual={:.3e} kg",
        inventory.dry_solids()?.value(),
        fed - discharged - inventory.dry_solids()?.value(),
    );
    Ok(())
}
