//! Published Co-60 -> Ni-60 with energy-resolved interval gamma emission.
//! Run: cargo run -p fs-nuclear --example evaluated_decay_emissions
//! Original evaluations and Git blobs: tests/data/endf_b_viii_0/README.md.
//!
//! The original inventory owns populations and calorimetry. The prepared
//! spectrum is a separate reusable emission observable, not additional heat.

use std::sync::Arc;

use fs_la::decay::DecayLimits;
use fs_nuclear::decay::{DecayBudget, DecayInventory};
use fs_nuclear::endf::decay::model::emissions::{DecayEmissionLimits, PreparedDecayEmission};
use fs_nuclear::endf::decay::model::{
    DecayEnergyFractions, DecayEnergyRouting, DecayModelEntry, EvaluatedDecayAssembly,
};
use fs_nuclear::endf::decay::spectra::DecayRadiation;
use fs_nuclear::endf::decay::{DecayEvaluation, EndfDecayReadLimits, EndfDecaySelection};
use fs_nuclear::irradiation::NuclideId;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cobalt = NuclideId::new(27, 60, 0)?;
    let nickel = NuclideId::new(28, 60, 0)?;
    let inputs = [
        (
            "ENDF/B-VIII.0/Co60/blob-f62ee6f48d10dc78c0f0eae8620bc0232cd4fbab",
            558,
            cobalt,
            include_str!("../tests/data/endf_b_viii_0/dec-027_Co_060.endf"),
        ),
        (
            "ENDF/B-VIII.0/Ni60/blob-9f67cefc53b9a2bd2da63ddad768c1d825dc82c0",
            590,
            nickel,
            include_str!("../tests/data/endf_b_viii_0/dec-028_Ni_060.endf"),
        ),
    ];
    // Declared local stopping assumption for charged particles, independent
    // of the photon spectrum. No Q remainder supplies missing neutrino data.
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
        include_reported_neutrinos: false,
    };
    let mut entries = Vec::new();
    for (data_id, material_number, nuclide, text) in inputs {
        entries.push(DecayModelEntry {
            evaluation: Arc::new(DecayEvaluation::parse(
                data_id.into(),
                EndfDecaySelection {
                    material_number,
                    nuclide,
                },
                text,
                EndfDecayReadLimits::default(),
                || Ok(()),
            )?),
            routing,
        });
    }
    let limits = DecayLimits::default();
    let assembly = Arc::new(EvaluatedDecayAssembly::new(
        "ENDF/B-VIII.0/cobalt60-emission-example".into(),
        entries,
        limits,
        || Ok(()),
    )?);
    let parent = assembly.index(cobalt).ok_or("missing cobalt-60")?;
    let daughter = assembly.index(nickel).ok_or("missing nickel-60")?;
    let half_life_s = assembly.evaluations()[parent]
        .law()
        .half_life()
        .ok_or("cobalt-60 must have a finite half-life")?
        .value();
    const EV_J: f64 = 1.602_176_634e-19;
    let source = PreparedDecayEmission::new(
        Arc::clone(&assembly),
        DecayRadiation::Gamma,
        vec![0.0, 1.25e6 * EV_J, 3.0e6 * EV_J],
        DecayEmissionLimits::default(),
        || Ok(()),
    )?;
    let initial_parent = 1e12;
    let mut populations = vec![0.0; assembly.nuclides().len()];
    populations[parent] = initial_parent;
    let mut inventory = DecayInventory::new(Arc::clone(assembly.model()), populations)?;
    let prepared = assembly
        .model()
        .prepare(half_life_s / 4.0, limits, || Ok(()))?;
    let budget = DecayBudget {
        propagation: limits,
        nuclei_tolerance: 0.1,
        energy_tolerance_j: 1e-12,
    };
    let mut total_photons = 0.0;
    println!(
        "half_lives,cobalt_nuclei,nickel_nuclei,interval_low_gamma_photons,interval_high_gamma_photons,interval_gamma_energy_j,interval_deposited_j"
    );
    for sample in 1..=8 {
        let step = inventory.advance_prepared(&prepared, budget, || Ok(()))?;
        let emitted = source.from_decay_step(inventory.model(), &step, || Ok(()))?;
        assert!(emitted.unreported_decays.is_empty());
        total_photons += emitted.emitted_particles.iter().sum::<f64>();
        println!(
            "{:.2},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e}",
            sample as f64 / 4.0,
            inventory.populations()[parent],
            inventory.populations()[daughter],
            emitted.emitted_particles[0],
            emitted.emitted_particles[1],
            emitted.emitted_energy_j.iter().sum::<f64>(),
            step.released.deposited_j
        );
    }
    // Independent checks: two half-lives leave one quarter of the parent;
    // the six published gamma intensities sum to 1.99848902 photons per decay.
    assert!((inventory.populations()[parent] / initial_parent - 0.25).abs() < 1e-13);
    assert!((total_photons / (0.75 * initial_parent) - 1.99848902).abs() < 1e-12);
    Ok(())
}
