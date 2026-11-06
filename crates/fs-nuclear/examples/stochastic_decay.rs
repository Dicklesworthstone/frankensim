//! Eight manufactured atoms: a 2 s parent, 0.5 s daughter and stable endpoint.
//! Energies are authored synthetic joules/event, not an evaluated isotope table.
//! Compare integer sampled histories with the deterministic `decay_chain` example.
//! Run with `cargo run -p fs-nuclear --example stochastic_decay --locked`.

use std::sync::Arc;

use fs_la::decay::DecayLimits;
use fs_nuclear::EnergyChannels;
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
    let energy = |deposited_j, radiation_j| EnergyChannels {
        deposited_j,
        radiation_j,
        escaped_j: 0.0,
    };
    let branches = vec![
        DecayBranch {
            parent: 0,
            daughter: 1,
            probability: 1.0,
            energy_per_decay: energy(2e-13, 1e-13),
        },
        DecayBranch {
            parent: 1,
            daughter: 2,
            probability: 1.0,
            energy_per_decay: energy(1e-13, 2e-13),
        },
    ];
    let model = Arc::new(DecayModel::new(
        "manufactured-half-life-history-v1".into(),
        nodes,
        branches,
        DecayLimits::default(),
    )?);
    let budget = StochasticDecayBudget {
        maximum_events: 32,
        energy_tolerance_j: 1e-24,
    };
    let key = StreamKey {
        seed: 20261105,
        kernel: 0xdec0_0001,
        tile: 0,
    };
    let mut history = StochasticDecayInventory::new(model, vec![8, 0, 0], key)?;
    let mut whole = history.clone();
    whole.advance_to(12.0, budget, || Ok(()))?;

    println!("time_s  [parent, daughter, stable]  events  deposited_interval_j");
    let mut checkpoint = None;
    for end in [0.0, 0.5, 1.0, 2.0, 2.5, 3.0, 4.0, 6.0, 8.0, 12.0] {
        let report = history.advance_to(end, budget, || Ok(()))?;
        println!(
            "{end:6.2}  {:?}  {:6}  {:.8e}",
            history.populations(),
            report.events,
            report.released.deposited_j,
        );
        // Retain the physical state, pending clocks and RNG together.
        if end == 2.0 {
            checkpoint = Some(history.clone());
        }
    }
    let mut restarted = checkpoint.expect("the 2 s checkpoint was retained");
    restarted.advance_to(12.0, budget, || Ok(()))?;
    assert_eq!(history, whole);
    assert_eq!(restarted, whole);
    println!("Partitioned observation and cloned restart retain the identical seeded trajectory.");
    println!(
        "Lifetime released energy [J]: {:?}",
        history.cumulative_released()
    );
    Ok(())
}
