//! Published ENDF/B-VIII.0 tritium -> stable helium-3, using the original decay owner.
//! Run: cargo run -p fs-nuclear --example evaluated_decay
//! Fixture source and identities: tests/data/endf_b_viii_0/README.md.
//!
//! The explicit example assumption deposits all electron/heavy-particle energy
//! and exports electromagnetic energy as radiation. The imported half-life and
//! mean emission energies are evaluation central values. Their uncertainties
//! remain available on the evaluation; they are not sampled in this example.

use std::sync::Arc;

use fs_la::decay::DecayLimits;
use fs_nuclear::decay::{DecayBudget, DecayInventory};
use fs_nuclear::endf::decay::model::{
    DecayEnergyFractions, DecayEnergyRouting, DecayModelEntry, EvaluatedDecayAssembly,
};
use fs_nuclear::endf::decay::{DecayEvaluation, EndfDecayReadLimits, EndfDecaySelection};
use fs_nuclear::irradiation::NuclideId;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let tritium = NuclideId::new(1, 3, 0)?;
    let helium = NuclideId::new(2, 3, 0)?;
    let inputs = [
        (
            "ENDF/B-VIII.0/H3/blob-306f60015644d4ba758217d6b4e66ac96986998f",
            EndfDecaySelection {
                material_number: 4,
                nuclide: tritium,
            },
            include_str!("../tests/data/endf_b_viii_0/dec-001_H_003.endf"),
        ),
        (
            "ENDF/B-VIII.0/He3/blob-6941a4d431ca4a9b48a77d8fad807a695c84f957",
            EndfDecaySelection {
                material_number: 9,
                nuclide: helium,
            },
            include_str!("../tests/data/endf_b_viii_0/dec-002_He_003.endf"),
        ),
    ];
    let routing = DecayEnergyRouting {
        light_particles: DecayEnergyFractions {
            deposited: 1.0,
            radiation: 0.0,
            escaped: 0.0,
        },
        electromagnetic: DecayEnergyFractions {
            deposited: 0.0,
            radiation: 1.0,
            escaped: 0.0,
        },
        heavy_particles: DecayEnergyFractions {
            deposited: 1.0,
            radiation: 0.0,
            escaped: 0.0,
        },
        // NC=3 has no reported neutrino category. Do not fill it from Q.
        include_reported_neutrinos: false,
    };
    let mut entries = Vec::new();
    for (id, selection, text) in inputs {
        let evaluation = DecayEvaluation::parse(
            id.into(),
            selection,
            text,
            EndfDecayReadLimits::default(),
            || Ok(()),
        )?;
        entries.push(DecayModelEntry {
            evaluation: Arc::new(evaluation),
            routing,
        });
    }
    let limits = DecayLimits::default();
    let assembly = EvaluatedDecayAssembly::new(
        "ENDF/B-VIII.0/tritium-to-helium3/example".into(),
        entries,
        limits,
        || Ok(()),
    )?;
    let parent = assembly.index(tritium).ok_or("missing tritium")?;
    let daughter = assembly.index(helium).ok_or("missing helium-3")?;
    let half_life_s = assembly.evaluations()[parent]
        .law()
        .half_life()
        .ok_or("tritium must have a finite half-life")?
        .value();
    let step_s = half_life_s / 4.0;
    let budget = DecayBudget {
        propagation: limits,
        nuclei_tolerance: 0.1,
        energy_tolerance_j: 1e-12,
    };
    let mut initial = vec![0.0; assembly.nuclides().len()];
    initial[parent] = 1e12;
    let mut inventory = DecayInventory::new(Arc::clone(assembly.model()), initial)?;
    let prepared = assembly.model().prepare(step_s, limits, || Ok(()))?;

    println!(
        "time_s,half_lives,tritium_nuclei,helium3_nuclei,tritium_activity_bq,cumulative_deposited_j,cumulative_radiation_j"
    );
    for sample in 0..=8 {
        if sample > 0 {
            inventory.advance_prepared(&prepared, budget, || Ok(()))?;
        }
        let populations = inventory.populations();
        let activities = inventory.activity_bq()?;
        let released = inventory.cumulative_released();
        println!(
            "{:.9e},{:.2},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e}",
            sample as f64 * step_s,
            sample as f64 / 4.0,
            populations[parent],
            populations[daughter],
            activities[parent],
            released.deposited_j,
            released.radiation_j,
        );
    }
    Ok(())
}
