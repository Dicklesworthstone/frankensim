//! A manufactured Poisson feed enters compartment A, moves to B, then drains.
//! Parent and daughter half-lives are 2 s and 0.5 s; the feed stops at 2 s.
//! Coefficients and energies are synthetic, not evaluated isotope or flow data.
//! Run with `cargo run -p fs-nuclear --example stochastic_production_transport --locked`.

use std::sync::Arc;

use fs_la::decay::DecayLimits;
use fs_nuclear::EnergyChannels;
use fs_nuclear::decay::production::{DecayProductionChange, DecayProductionSchedule};
use fs_nuclear::decay::stochastic::transport::production::{
    StochasticDecayProductionTransportBudget, StochasticDecayProductionTransportState,
};
use fs_nuclear::decay::stochastic::transport::{
    NuclideOutflow, NuclideTransfer, NuclideTransport, StochasticDecayTransportBudget,
    StochasticDecayTransportLimits, StochasticDecayTransportState,
};
use fs_nuclear::decay::stochastic::{StochasticDecayBudget, StochasticDecayInventory, StreamKey};
use fs_nuclear::decay::{DecayBranch, DecayModel, DecayNuclide};
use fs_qty::Time;
use fs_thermochem::decay::DecayLawV1;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let model = Arc::new(DecayModel::new(
        "manufactured-fed-compartment-half-life-history-v1".into(),
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
                    kernel: 0xdec0_0003,
                    tile: cell,
                },
            )
        })
        .collect::<Result<Vec<_>, _>>()?;
    let transport = StochasticDecayTransportState::new(
        inventories,
        NuclideTransport {
            transfers: vec![NuclideTransfer {
                from_compartment: 0,
                to_compartment: 1,
                rates_per_s: vec![1.5; 3],
            }],
            outflows: vec![NuclideOutflow {
                compartment: 1,
                rates_per_s: vec![0.15; 3],
            }],
        },
        StreamKey {
            seed: 20261106,
            kernel: 0xdec0_1003,
            tile: 0,
        },
        StochasticDecayTransportLimits::default(),
    )?;
    let schedule = Arc::new(DecayProductionSchedule::new(
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
        StochasticDecayProductionTransportState::poisson(transport, vec![Some(schedule), None])?;
    let mut whole = history.clone();
    whole.advance_to(20.0, budget, || Ok(()))?;

    println!("time_s  A[P,D,S]  B[P,D,S]  born  decays  outflow_total  heat_A_J  heat_B_J");
    let mut checkpoint = None;
    for end in [0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 20.0] {
        let step = history.advance_to(end, budget, || Ok(()))?;
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
            .expect("A has the explicit source")
            .iter()
            .sum::<u64>();
        // One represented daughter per decay; internal transfers do not create nuclei.
        assert_eq!(supplied, retained + outflow);
        println!(
            "{end:6.2}  {:?}  {:?}  {:4}  {:6}  {:13}  {:.6e}  {:.6e}",
            history.inventories()[0].populations(),
            history.inventories()[1].populations(),
            step.births,
            step.transport.decay_events,
            outflow,
            step.transport.decay[0].released.deposited_j,
            step.transport.decay[1].released.deposited_j,
        );
        if end == 1.0 {
            // Retain the source and movement clocks together with both inventories.
            checkpoint = Some(history.clone());
        }
    }
    let mut restarted = checkpoint.expect("the 1 s checkpoint was retained");
    restarted.advance_to(20.0, budget, || Ok(()))?;
    assert_eq!(history, whole);
    assert_eq!(restarted, whole);
    println!("Exact integer balance holds at every observation: supplied = retained + outflow.");
    println!(
        "Births and transfers add no inferred heat; emissions remain at their event compartment."
    );
    println!(
        "Partitioned observation and cloned restart retain the identical complete trajectory."
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
