//! Published Co-60 -> Ni-60 with energy-resolved gamma emission and shielding.
//! Run: cargo run -p fs-nuclear --example evaluated_decay_emissions
//! Original evaluations and Git blobs: tests/data/endf_b_viii_0/README.md.
//!
//! The original inventory owns populations and calorimetry. A manufactured
//! planar photon kernel demonstrates multiple scattering and lower-group
//! energy transfer. These illustrative coefficients are NOT evaluated material
//! shielding data. Deposited photon energy is drawn from the original outward
//! radiation account; the spectral observable creates no additional energy.

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
use fs_nuclear::radiation::decay::{DecayRadiationClosure, PreparedDecayRadiation};
use fs_nuclear::radiation::transport::refinement::{PhotonSlabRefinement, refine_slab_response};
use fs_nuclear::radiation::transport::{
    PhotonCollision, PhotonSlabCell, PhotonSlabLimits, gauss_legendre_4,
};

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
    let source = Arc::new(PreparedDecayEmission::new(
        Arc::clone(&assembly),
        DecayRadiation::Gamma,
        vec![0.0, 1.25e6 * EV_J, 3.0e6 * EV_J],
        DecayEmissionLimits::default(),
        || Ok(()),
    )?);
    // An explicitly declared allowance for the evaluated spectrum/aggregate
    // mean-energy discrepancy, not a transport or numerical-error allowance.
    let radiative = PreparedDecayRadiation::new(
        Arc::clone(&source),
        DecayRadiationClosure {
            absolute_j_per_decay: 0.0,
            relative: 1e-6,
        },
        DecayEmissionLimits::default(),
        || Ok(()),
    )?;
    let shielding = refine_slab_response(
        "manufactured two-group isotropic planar shielding example".into(),
        DecayRadiation::Gamma,
        source.group_edges_j().to_vec(),
        vec![
            PhotonSlabCell {
                cell_id: 10,
                thickness_m: 0.001,
                // Thin-source transparent approximation, declared explicitly.
                groups: vec![
                    PhotonCollision {
                        total_per_m: 0.0,
                        heat_per_m: 0.0,
                        redistribution_per_m: vec![0.0, 0.0],
                    };
                    2
                ],
            },
            PhotonSlabCell {
                cell_id: 20,
                thickness_m: 0.02,
                groups: vec![
                    PhotonCollision {
                        total_per_m: 50.0,
                        heat_per_m: 35.0,
                        redistribution_per_m: vec![15.0, 0.0],
                    },
                    PhotonCollision {
                        total_per_m: 30.0,
                        heat_per_m: 5.0,
                        redistribution_per_m: vec![15.0, 10.0],
                    },
                ],
            },
        ],
        vec![10],
        gauss_legendre_4(),
        PhotonSlabLimits::default(),
        PhotonSlabRefinement {
            absolute_tolerance: 1e-4,
            relative_tolerance: 1e-3,
            max_levels: 6,
        },
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
        "half_lives,cobalt_nuclei,nickel_nuclei,interval_low_gamma_photons,interval_high_gamma_photons,interval_gamma_energy_j,interval_deposited_j,interval_shield_heat_j,interval_photon_escape_j,interval_routed_gamma_energy_j,interval_spectral_adjustment_j"
    );
    for sample in 1..=8 {
        let step = inventory.advance_prepared(&prepared, budget, || Ok(()))?;
        let emitted = source.from_decay_step(inventory.model(), &step, || Ok(()))?;
        let partition = radiative.from_decay_step(inventory.model(), &step, || Ok(()))?;
        let absorbed = shielding.transport.response.apply(
            &[(10, &partition.emitted_energy_j)],
            1e-12,
            || Ok(()),
        )?;
        assert!(
            (absorbed.deposited_j + absorbed.escaped_j - step.released.radiation_j).abs() < 1e-12
        );
        assert!(emitted.unreported_decays.is_empty());
        total_photons += emitted.emitted_particles.iter().sum::<f64>();
        // The evaluated spectral mean and reported mean energy can differ.
        // Expose that admitted adjustment separately; shield plus escape closes
        // the original routed energy, and creates no additional decay heat.
        println!(
            "{:.2},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e}",
            sample as f64 / 4.0,
            inventory.populations()[parent],
            inventory.populations()[daughter],
            emitted.emitted_particles[0],
            emitted.emitted_particles[1],
            emitted.emitted_energy_j.iter().sum::<f64>(),
            step.released.deposited_j,
            absorbed.deposited_j,
            absorbed.escaped_j,
            partition.original_radiation_j,
            partition.reconciled_spectral_energy_j
        );
    }
    // Independent checks: two half-lives leave one quarter of the parent;
    // the six published gamma intensities sum to 1.99848902 photons per decay.
    assert!((inventory.populations()[parent] / initial_parent - 0.25).abs() < 1e-13);
    assert!((total_photons / (0.75 * initial_parent) - 1.99848902).abs() < 1e-12);
    Ok(())
}
