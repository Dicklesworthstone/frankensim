//! Synthetic production, daughter ingrowth and shutdown using one prepared map.
//! Run: cargo run -p fs-nuclear --example decay_production
//! All nuclear data and production rates below are authored demonstration values.

use fs_la::decay::DecayLimits;
use fs_nuclear::{
    EnergyChannels,
    decay::{
        DecayBranch, DecayBudget, DecayInventory, DecayModel, DecayNuclide,
        production::{DecayProductionChange, DecayProductionSchedule},
    },
};
use fs_qty::Time;
use fs_thermochem::decay::DecayLawV1;
use std::sync::Arc;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let nuclides = [
        ("synthetic-parent", Some(2.0)),
        ("synthetic-daughter", Some(0.5)),
        ("synthetic-stable", None),
    ]
    .into_iter()
    .map(|(id, half_life)| {
        Ok(DecayNuclide {
            id: id.into(),
            law: match half_life {
                Some(seconds) => DecayLawV1::from_half_life(Time::new(seconds))?,
                None => DecayLawV1::stable(),
            },
        })
    })
    .collect::<Result<Vec<_>, Box<dyn std::error::Error>>>()?;
    let branches = vec![
        DecayBranch {
            parent: 0,
            daughter: 1,
            probability: 1.0,
            energy_per_decay: EnergyChannels {
                deposited_j: 2e-13,
                radiation_j: 1e-13,
                escaped_j: 0.5e-13,
            },
        },
        DecayBranch {
            parent: 1,
            daughter: 2,
            probability: 1.0,
            energy_per_decay: EnergyChannels {
                deposited_j: 3e-13,
                radiation_j: 0.25e-13,
                escaped_j: 0.75e-13,
            },
        },
    ];
    let limits = DecayLimits::default();
    let model = Arc::new(DecayModel::new(
        "synthetic-production-shutdown/v1".into(),
        nuclides,
        branches,
        limits,
    )?);
    let schedule = DecayProductionSchedule::new(
        Arc::clone(&model),
        vec![1e12, 0.0, 0.0],
        vec![
            DecayProductionChange {
                time_s: 0.75,
                rates_nuclei_s: vec![2e12, 0.0, 0.0],
            },
            DecayProductionChange {
                time_s: 1.5,
                rates_nuclei_s: vec![0.0; 3],
            },
        ],
        2,
    )?;
    let budget = DecayBudget {
        propagation: limits,
        nuclei_tolerance: 1.0,
        energy_tolerance_j: 1e-12,
    };
    let dt = 0.25;
    let prepared = model.prepare_production(dt, limits, || Ok(()))?;
    let mut inventory = DecayInventory::new(model, vec![0.0; 3])?;
    let mut shutdown_heat = 0.0;
    println!(
        "time_s,parent_nuclei,daughter_nuclei,stable_nuclei,produced_nuclei,deposited_j,radiation_j,escaped_j"
    );
    for step in 0..16 {
        let start_s = step as f64 * dt;
        // This example's transitions fall exactly on the chosen grid. The
        // thermal specialization clips arbitrary authored transitions itself.
        assert!(
            schedule
                .next_after(start_s)
                .is_none_or(|next| next >= start_s + dt)
        );
        let report = inventory.advance_production_prepared(
            &prepared,
            schedule.rates_at(start_s)?,
            budget,
            || Ok(()),
        )?;
        let produced = inventory.cumulative_produced_nuclei().iter().sum::<f64>();
        let populations = inventory.populations();
        assert!((populations.iter().sum::<f64>() - produced).abs() < 2.0);
        let released = inventory.cumulative_released();
        if start_s >= 1.5 {
            assert!(report.produced_nuclei.iter().all(|&value| value == 0.0));
            shutdown_heat += report.decay.released.deposited_j;
        }
        println!(
            "{:.2},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e}",
            start_s + dt,
            populations[0],
            populations[1],
            populations[2],
            produced,
            released.deposited_j,
            released.radiation_j,
            released.escaped_j
        );
    }
    assert!((inventory.cumulative_produced_nuclei()[0] - 2.25e12).abs() < 1.0);
    assert!(shutdown_heat > 0.0);
    assert!(inventory.cumulative_decays()[1] > 0.0);
    Ok(())
}
