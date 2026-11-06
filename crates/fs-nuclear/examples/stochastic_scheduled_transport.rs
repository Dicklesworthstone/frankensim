//! Manufactured Poisson feed, valve changes, extraction and daughter half-lives.
//! A and B are fixed well-mixed 0.01 and 0.02 m3 control volumes. Each path uses
//! its explicit volumetric flow divided by its donor volume. These coefficients
//! and energies are illustrative, not evaluated isotope or calibrated flow data.
//! Run with `cargo run -p fs-nuclear --example stochastic_scheduled_transport --locked`.

use fs_la::decay::DecayLimits;
use fs_nuclear::burnup::transport::schedule::NuclideTransportChange;
use fs_nuclear::decay::production::{DecayProductionChange, DecayProductionSchedule};
use fs_nuclear::decay::stochastic::transport::production::{
    StochasticDecayProductionTransportBudget, StochasticDecayProductionTransportState,
};
use fs_nuclear::decay::stochastic::transport::schedule::StochasticDecayTransportSchedule;
use fs_nuclear::decay::stochastic::transport::{
    NuclideOutflow, NuclideTransfer, NuclideTransport, StochasticDecayTransportBudget,
    StochasticDecayTransportLimits, StochasticDecayTransportState,
};
use fs_nuclear::decay::stochastic::{StochasticDecayBudget, StochasticDecayInventory, StreamKey};
use fs_nuclear::decay::{DecayBranch, DecayModel, DecayNuclide};
use fs_nuclear::{EnergyChannels, FissionError};
use fs_qty::{Time, Volume, VolumetricFlowRate};
use fs_thermochem::decay::DecayLawV1;
use std::sync::Arc;

fn frame(
    volumes: &[Volume],
    transfer_m3_s: f64,
    extraction_m3_s: f64,
) -> Result<NuclideTransport, FissionError> {
    Ok(NuclideTransport {
        transfers: vec![NuclideTransfer::from_volumetric_flow(
            0,
            1,
            volumes,
            VolumetricFlowRate::new(transfer_m3_s),
            3,
        )?],
        outflows: vec![NuclideOutflow::from_volumetric_flow(
            1,
            volumes,
            VolumetricFlowRate::new(extraction_m3_s),
            3,
        )?],
    })
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let model = Arc::new(DecayModel::new(
        "manufactured-scheduled-flow-half-life-history-v1".into(),
        vec![
            DecayNuclide {
                id: "manufactured-parent".into(),
                law: DecayLawV1::from_half_life(Time::new(2.0))?,
            },
            DecayNuclide {
                id: "manufactured-daughter".into(),
                law: DecayLawV1::from_half_life(Time::new(0.5))?,
            },
            DecayNuclide {
                id: "manufactured-stable".into(),
                law: DecayLawV1::stable(),
            },
        ],
        vec![
            DecayBranch {
                parent: 0,
                daughter: 1,
                probability: 1.0,
                energy_per_decay: EnergyChannels {
                    deposited_j: 2e-13,
                    radiation_j: 1e-13,
                    escaped_j: 0.0,
                },
            },
            DecayBranch {
                parent: 1,
                daughter: 2,
                probability: 1.0,
                energy_per_decay: EnergyChannels {
                    deposited_j: 1e-13,
                    radiation_j: 2e-13,
                    escaped_j: 0.0,
                },
            },
        ],
        DecayLimits::default(),
    )?);
    let inventories = (0..2)
        .map(|cell| {
            StochasticDecayInventory::new(
                Arc::clone(&model),
                vec![0; 3],
                StreamKey {
                    seed: 20261106,
                    kernel: 0xdec0_0004,
                    tile: cell,
                },
            )
        })
        .collect::<Result<Vec<_>, _>>()?;
    let volumes = [Volume::new(0.01), Volume::new(0.02)];
    let closed = frame(&volumes, 0.0, 0.0)?;
    let transport = StochasticDecayTransportState::new(
        inventories,
        closed.clone(),
        StreamKey {
            seed: 20261106,
            kernel: 0xdec0_1004,
            tile: 0,
        },
        StochasticDecayTransportLimits::default(),
    )?;
    let movement = Arc::new(StochasticDecayTransportSchedule::new(
        &transport,
        closed,
        vec![
            NuclideTransportChange {
                time_s: 0.5,
                transport: frame(&volumes, 0.015, 0.0)?,
            },
            NuclideTransportChange {
                time_s: 2.0,
                transport: frame(&volumes, 0.0, 0.003)?,
            },
            NuclideTransportChange {
                time_s: 3.0,
                transport: frame(&volumes, 0.01, 0.003)?,
            },
        ],
        3,
    )?);
    let supply = Arc::new(DecayProductionSchedule::new(
        model,
        vec![4.0, 0.0, 0.0],
        vec![DecayProductionChange {
            time_s: 2.0,
            rates_nuclei_s: vec![0.0; 3],
        }],
        1,
    )?);
    let budget = StochasticDecayProductionTransportBudget {
        transport: StochasticDecayTransportBudget {
            decay: StochasticDecayBudget {
                maximum_events: 256,
                energy_tolerance_j: 1e-24,
            },
            maximum_movements: 512,
            maximum_clock_checks: 100_000,
        },
        maximum_births: 128,
        maximum_schedule_changes: 1,
    };
    let mut history =
        StochasticDecayProductionTransportState::poisson(transport, vec![Some(supply), None])?
            .with_transport_schedule(movement)?;
    let mut whole = history.clone();
    whole.advance_to(20.0, budget, || Ok(()))?;

    println!("Feed: 4 parents/s into A until 2 s; parent/daughter half-lives: 2/0.5 s.");
    println!("Valve opens at 0.5 s, closes at 2 s as B extraction begins, and reopens at 3 s.");
    println!(
        "time_s  A[P,D,S]  B[P,D,S]  born  decays  flow_changes  outflow_total  heat_A_J  heat_B_J"
    );
    let mut checkpoint = None;
    let mut changes = 0;
    for end in [
        0.0, 0.25, 0.5, 1.0, 2.0, 2.5, 3.0, 4.0, 6.0, 8.0, 12.0, 20.0,
    ] {
        let step = history.advance_to(end, budget, || Ok(()))?;
        changes += step.transport.transport_changes;
        let outflow = history
            .transport_state()
            .cumulative_outflow_counts()
            .iter()
            .flatten()
            .sum::<u64>();
        let retained = history
            .inventories()
            .iter()
            .flat_map(|cell| cell.populations())
            .sum::<u64>();
        let supplied = history
            .cumulative_produced_nuclei(0)
            .unwrap()
            .iter()
            .sum::<u64>();
        assert_eq!(supplied, retained + outflow);
        if end > 2.0 {
            assert_eq!(step.births, 0);
        }
        println!(
            "{end:6.2}  {:?}  {:?}  {:4}  {:6}  {:12}  {:13}  {:.6e}  {:.6e}",
            history.inventories()[0].populations(),
            history.inventories()[1].populations(),
            step.births,
            step.transport.decay_events,
            step.transport.transport_changes,
            outflow,
            step.transport.decay[0].released.deposited_j,
            step.transport.decay[1].released.deposited_j
        );
        if end == 1.0 {
            checkpoint = Some(history.clone());
        }
    }
    let mut restarted = checkpoint.expect("the open-valve, active-source checkpoint was retained");
    restarted.advance_to(20.0, budget, || Ok(()))?;
    assert_eq!(history, whole);
    assert_eq!(restarted, whole);
    assert_eq!(changes, 3);
    println!(
        "Supplied = retained + outflow at every observation; all three valve changes apply once."
    );
    println!(
        "Partitioned observation and cloned restart preserve sources, valves, clocks and events."
    );
    println!(
        "Births and carrier movement add no inferred heat; only local decay deposits its branch energy."
    );
    println!(
        "Lifetime supplied [P,D,S]: {:?}",
        history.cumulative_produced_nuclei(0).unwrap()
    );
    println!(
        "Lifetime outflow [P,D,S]: {:?}",
        history.transport_state().cumulative_outflow_counts()[0]
    );
    for (cell, inventory) in history.inventories().iter().enumerate() {
        println!(
            "Compartment {cell} lifetime released [J]: {:?}",
            inventory.cumulative_released()
        );
    }
    Ok(())
}
