//! Finite targets, repeated capture, daughter half-lives and a shutdown tail.
//!
//! All coefficients below are manufactured equation data. The nuclide identities
//! make the reaction products explicit; the supplied seconds-long half-lives,
//! cross sections and energies are NOT evaluated data for these isotopes.
//! Run with `cargo run -p fs-nuclear --example stochastic_activation --locked`.

use std::sync::Arc;

use fs_nuclear::EnergyChannels;
use fs_nuclear::burnup::stochastic::{
    StochasticActivationBudget, StochasticActivationChange, StochasticActivationSchedule,
    StochasticActivationState, StreamKey,
};
use fs_nuclear::burnup::{
    BurnupData, BurnupLimits, BurnupModel, BurnupNuclide, BurnupProduct, BurnupRate, BurnupReaction,
};
use fs_nuclear::irradiation::{GroupFlux, NeutronGroups, NuclideId};
use fs_qty::Time;
use fs_thermochem::decay::DecayLawV1;

fn reaction(
    label: &str,
    parent: usize,
    daughter: usize,
    rate: BurnupRate,
    deposited_j: f64,
    radiation_j: f64,
) -> BurnupReaction {
    BurnupReaction {
        data_id: label.into(),
        parent,
        rate,
        products: vec![BurnupProduct {
            nuclide: daughter,
            yield_per_event: 1.0,
        }],
        energy_per_event: EnergyChannels {
            deposited_j,
            radiation_j,
            escaped_j: 0.0,
        },
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let groups = Arc::new(NeutronGroups::new(
        "manufactured two-group activation spectrum".into(),
        vec![0.0, 1e-20, 1e-13],
    )?);
    let stable = DecayLawV1::stable();
    let model = Arc::new(BurnupModel::new(
        BurnupData {
            data_id: "manufactured stochastic capture and daughter-recapture example".into(),
            groups: Arc::clone(&groups),
            temperature_min_k: 200.0,
            temperature_max_k: 1000.0,
            nuclides: vec![
                BurnupNuclide {
                    id: NuclideId::new(27, 59, 0)?,
                    law: stable,
                },
                BurnupNuclide {
                    id: NuclideId::new(27, 60, 0)?,
                    law: DecayLawV1::from_half_life(Time::new(2.0))?,
                },
                BurnupNuclide {
                    id: NuclideId::new(27, 61, 0)?,
                    law: DecayLawV1::from_half_life(Time::new(0.5))?,
                },
                BurnupNuclide {
                    id: NuclideId::new(28, 60, 0)?,
                    law: stable,
                },
                BurnupNuclide {
                    id: NuclideId::new(28, 61, 0)?,
                    law: stable,
                },
            ],
            reactions: vec![
                reaction(
                    "target capture",
                    0,
                    1,
                    BurnupRate::Capture {
                        cross_sections_m2: vec![1e-28, 2e-29],
                    },
                    2e-13,
                    1e-13,
                ),
                reaction(
                    "activated daughter capture",
                    1,
                    2,
                    BurnupRate::Capture {
                        cross_sections_m2: vec![1e-29, 4e-29],
                    },
                    1e-13,
                    2e-13,
                ),
                reaction(
                    "first daughter decay",
                    1,
                    3,
                    BurnupRate::Decay {
                        branch_fraction: 1.0,
                    },
                    5e-14,
                    2e-14,
                ),
                reaction(
                    "second daughter decay",
                    2,
                    4,
                    BurnupRate::Decay {
                        branch_fraction: 1.0,
                    },
                    8e-14,
                    1e-14,
                ),
            ],
        },
        BurnupLimits::default(),
    )?);
    let flux = |values| GroupFlux::new(Arc::clone(&groups), values);
    let schedule = Arc::new(StochasticActivationSchedule::new(
        model,
        flux(vec![3e28, 1e27])?,
        300.0,
        vec![
            StochasticActivationChange {
                time_s: 1.5,
                flux: flux(vec![0.5e28, 2e28])?,
                temperature_k: 300.0,
            },
            StochasticActivationChange {
                time_s: 3.0,
                flux: flux(vec![0.0; 2])?,
                temperature_k: 300.0,
            },
        ],
        BurnupLimits::default(),
        2,
        || Ok(()),
    )?);
    let mut history = StochasticActivationState::new(
        schedule,
        vec![64, 0, 0, 0, 0],
        StreamKey {
            seed: 20261107,
            kernel: 0x6163_7469,
            tile: 1,
        },
    )?;
    let budget = StochasticActivationBudget {
        maximum_events: 1024,
        maximum_clock_checks: 100_000,
        energy_tolerance_j: 1e-24,
    };
    let mut whole = history.clone();
    whole.advance_to(100.0, budget, || Ok(()))?;
    println!("Manufactured five-nuclide network; daughter half-lives 2 s and 0.5 s.");
    println!("Spectrum changes at 1.5 s; irradiation stops at 3 s.");
    println!(
        "time_s populations[target,daughter1,daughter2,stable1,stable2] captures decays capture_heat_J decay_heat_J"
    );
    let mut shutdown = None;
    for end_s in [0.0, 0.25, 0.75, 1.5, 3.0, 4.0, 6.0, 10.0, 20.0, 100.0] {
        history.advance_to(end_s, budget, || Ok(()))?;
        assert_eq!(history.populations().iter().sum::<u64>(), 64);
        let counts = history.cumulative_reaction_counts();
        println!(
            "{end_s:6.2} {:?} {} {} {:.9e} {:.9e}",
            history.populations(),
            counts[0] + counts[1],
            counts[2] + counts[3],
            history.cumulative_capture_released().deposited_j,
            history.cumulative_decay_released().deposited_j,
        );
        if end_s == 3.0 {
            shutdown = Some(history.clone());
        }
    }
    assert_eq!(
        history, whole,
        "observations must preserve the physical trajectory"
    );
    let mut restarted = shutdown.expect("explicit shutdown checkpoint");
    let shutdown_captures = restarted.cumulative_reaction_counts()[..2].to_vec();
    let shutdown_decay_heat = restarted.cumulative_decay_released().deposited_j;
    restarted.advance_to(100.0, budget, || Ok(()))?;
    assert_eq!(
        restarted, history,
        "restart must preserve pending clocks and RNG"
    );
    assert_eq!(
        &history.cumulative_reaction_counts()[..2],
        &shutdown_captures
    );
    assert!(history.cumulative_decay_released().deposited_j > shutdown_decay_heat);
    assert_eq!(&history.populations()[1..3], &[0, 0]);
    println!("Exact 64-nucleus conservation and checkpoint/restart agreement verified.");
    Ok(())
}
