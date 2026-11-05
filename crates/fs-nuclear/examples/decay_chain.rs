//! Synthetic parent -> daughter -> stable inventory with integrated emission.
//! Run: cargo run -p fs-nuclear --example decay_chain
//! Half-lives and branch energies are explicit demonstration data, not an
//! empirical isotope model. Conserved tracked nuclei do not imply nuclear
//! mass/charge conservation without an emitted-particle inventory.

use fs_la::decay::DecayLimits;
use fs_nuclear::EnergyChannels;
use fs_nuclear::decay::{DecayBranch, DecayBudget, DecayInventory, DecayModel, DecayNuclide};
use fs_qty::Time;
use fs_thermochem::decay::DecayLawV1;
use std::sync::Arc;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let nuclides = vec![
        DecayNuclide {
            id: "synthetic-parent".into(),
            law: DecayLawV1::from_half_life(Time::new(2.0))?,
        },
        DecayNuclide {
            id: "synthetic-daughter".into(),
            law: DecayLawV1::from_half_life(Time::new(0.5))?,
        },
        DecayNuclide {
            id: "synthetic-stable".into(),
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
        "synthetic-three-node-chain/v1".into(),
        nuclides,
        branches,
        limits,
    )?);
    let budget = DecayBudget {
        propagation: limits,
        nuclei_tolerance: 0.1,
        energy_tolerance_j: 1e-12,
    };
    let dt = Time::new(0.5);
    let prepared = model.prepare(dt.value(), limits, || Ok(()))?;
    let initial_nuclei = 1e12;
    let mut inventory = DecayInventory::new(model, vec![initial_nuclei, 0.0, 0.0])?;

    println!(
        "time_s,parent_nuclei,daughter_nuclei,stable_nuclei,parent_activity_bq,daughter_activity_bq,stable_activity_bq,cumulative_deposited_j,cumulative_radiation_j,cumulative_escaped_j,cumulative_daughter_heat_j"
    );
    for step in 1..=16 {
        let report = inventory.advance_prepared(&prepared, budget, || Ok(()))?;
        let nuclei = inventory.populations();
        assert!((nuclei.iter().sum::<f64>() - initial_nuclei).abs() < 0.5);
        assert_eq!(report.branch_decays, report.decays[..2]);
        let activities = inventory.activity_bq()?;
        let released = inventory.cumulative_released();
        let daughter_heat_j = inventory.cumulative_decays()[1] * 3e-13;
        println!(
            "{:.1},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e}",
            step as f64 * dt.value(),
            nuclei[0],
            nuclei[1],
            nuclei[2],
            activities[0],
            activities[1],
            activities[2],
            released.deposited_j,
            released.radiation_j,
            released.escaped_j,
            daughter_heat_j
        );
    }
    let counts = inventory.cumulative_decays();
    let released = inventory.cumulative_released();
    assert!((inventory.populations()[0] - initial_nuclei / 16.0).abs() < 0.1);
    assert!((counts[0] - (initial_nuclei - inventory.populations()[0])).abs() < 0.5);
    assert!((counts[1] - inventory.populations()[2]).abs() < 0.5);
    assert!(counts[1] * 3e-13 > 0.0 && released.deposited_j > counts[1] * 3e-13);
    for (actual, expected) in [
        (released.deposited_j, counts[0] * 2e-13 + counts[1] * 3e-13),
        (
            released.radiation_j,
            counts[0] * 1e-13 + counts[1] * 0.25e-13,
        ),
        (
            released.escaped_j,
            counts[0] * 0.5e-13 + counts[1] * 0.75e-13,
        ),
    ] {
        assert!((actual - expected).abs() < 1e-12);
    }
    Ok(())
}
