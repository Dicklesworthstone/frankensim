//! Manufactured Poisson production, source shutdown, and two-stage radioactive decay.
//! Parent production is 4 nuclei/s until 2 s, then 2 nuclei/s from 4 to 6 s.
//! Half-lives and joules/event are synthetic, not an evaluated isotope table.
//! Run with `cargo run -p fs-nuclear --example stochastic_production --locked`.

use std::sync::Arc;

use fs_la::decay::DecayLimits;
use fs_nuclear::EnergyChannels;
use fs_nuclear::decay::production::{DecayProductionChange, DecayProductionSchedule};
use fs_nuclear::decay::stochastic::production::{
    StochasticDecayProduction, StochasticDecayProductionBudget,
};
use fs_nuclear::decay::stochastic::{StochasticDecayBudget, StochasticDecayInventory, StreamKey};
use fs_nuclear::decay::{DecayBranch, DecayModel, DecayNuclide};
use fs_qty::Time;
use fs_thermochem::decay::DecayLawV1;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let nodes = vec![
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
    ];
    let branches = vec![
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
    ];
    let model = Arc::new(DecayModel::new(
        "manufactured-poisson-half-life-history-v1".into(),
        nodes,
        branches,
        DecayLimits::default(),
    )?);
    let schedule = Arc::new(DecayProductionSchedule::new(
        Arc::clone(&model),
        vec![4.0, 0.0, 0.0],
        vec![
            DecayProductionChange {
                time_s: 2.0,
                rates_nuclei_s: vec![0.0; 3],
            },
            DecayProductionChange {
                time_s: 4.0,
                rates_nuclei_s: vec![2.0, 0.0, 0.0],
            },
            DecayProductionChange {
                time_s: 6.0,
                rates_nuclei_s: vec![0.0; 3],
            },
        ],
        3,
    )?);
    let inventory = StochasticDecayInventory::new(
        model,
        vec![0; 3],
        StreamKey {
            seed: 20261106,
            kernel: 0xdec0_0002,
            tile: 0,
        },
    )?;
    let budget = StochasticDecayProductionBudget {
        decay: StochasticDecayBudget {
            maximum_events: 256,
            energy_tolerance_j: 1e-24,
        },
        maximum_births: 128,
        maximum_schedule_changes: 3,
    };
    let mut history = StochasticDecayProduction::poisson(inventory, schedule)?;
    let mut whole = history.clone();
    whole.advance_to(20.0, budget, || Ok(()))?;

    println!("time_s  [parent, daughter, stable]  born  decays  deposited_interval_j");
    let mut checkpoint = None;
    for end in [0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 12.0, 20.0] {
        let step = history.advance_to(end, budget, || Ok(()))?;
        println!(
            "{end:6.2}  {:?}  {:4}  {:6}  {:.8e}",
            history.inventory().populations(),
            step.births,
            step.decay.events,
            step.decay.released.deposited_j,
        );
        // A restart retains pending births, decay clocks, schedule position and RNG.
        if end == 2.0 {
            checkpoint = Some(history.clone());
        }
    }
    let mut restarted = checkpoint.expect("the 2 s checkpoint was retained");
    restarted.advance_to(20.0, budget, || Ok(()))?;
    assert_eq!(history, whole);
    assert_eq!(restarted, whole);
    assert_eq!(
        history.inventory().populations().iter().sum::<u64>(),
        history.cumulative_produced_nuclei().iter().sum::<u64>(),
    );
    println!(
        "Births carry no inferred heat; the original branch releases energy at each sampled decay."
    );
    println!("Partitioned observation and cloned restart retain the identical seeded trajectory.");
    println!(
        "Lifetime produced nuclei: {:?}",
        history.cumulative_produced_nuclei()
    );
    println!(
        "Lifetime released energy [J]: {:?}",
        history.inventory().cumulative_released()
    );
    Ok(())
}
