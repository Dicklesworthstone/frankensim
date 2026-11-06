//! Evaluated tritium decay during prescribed circulation between two volumes.
//! Run: cargo run -p fs-nuclear --example evaluated_decay_transport
//!
//! ENDF/B-VIII.0 supplies the H3/He3 identities, half-life and emitted energy.
//! The circulation is a manufactured, well-mixed trace-inventory example:
//! V0 = 1 m³, V1 = 2 m³ and the same volumetric flow in both directions.
//! Q/V0 equals the imported decay constant so both processes are visible on
//! one half-life scale. These are authored flow assumptions, not plant data.
//! Every isotope has that donor-volume transfer rate. Carrier mass, momentum,
//! enthalpy and radiation transport remain outside this prescribed example.
//!
//! Stored reaction histories stay in their original compartments. The output
//! and independent analytic check distinguish current isotope location from
//! where earlier energy was released. Electron/heavy-particle energy is
//! explicitly deposited locally; electromagnetic energy leaves as radiation.

use std::sync::Arc;

use fs_la::decay::DecayLimits;
use fs_math::decay::{attenuate, decayed_fraction};
use fs_nuclear::burnup::transport::{
    BurnupTransportBudget, BurnupTransportLimits, BurnupTransportState, NuclideTransfer,
    NuclideTransport,
};
use fs_nuclear::burnup::{BurnupInventory, BurnupModel};
use fs_nuclear::endf::decay::model::{
    DecayEnergyFractions, DecayEnergyRouting, DecayModelEntry, EvaluatedDecayAssembly,
};
use fs_nuclear::endf::decay::{DecayEvaluation, EndfDecayReadLimits, EndfDecaySelection};
use fs_nuclear::irradiation::{GroupFlux, NeutronGroups, NuclideId};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let tritium = NuclideId::new(1, 3, 0)?;
    let helium = NuclideId::new(2, 3, 0)?;
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
        // The NC3 records provide no neutrino mean. Do not infer it from Q.
        include_reported_neutrinos: false,
    };
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
    let mut entries = Vec::new();
    for (data_id, selection, text) in inputs {
        entries.push(DecayModelEntry {
            evaluation: Arc::new(DecayEvaluation::parse(
                data_id.into(),
                selection,
                text,
                EndfDecayReadLimits::default(),
                || Ok(()),
            )?),
            routing,
        });
    }
    let assembly = EvaluatedDecayAssembly::new(
        "ENDF/B-VIII.0/H3-He3/prescribed-circulation-example".into(),
        entries,
        DecayLimits::default(),
        || Ok(()),
    )?;
    let parent = assembly.index(tritium).ok_or("missing tritium")?;
    let daughter = assembly.index(helium).ok_or("missing helium-3")?;
    let law = assembly.evaluations()[parent].law();
    let half_life_s = law.half_life().ok_or("finite half-life required")?.value();
    let lambda = law.decay_constant().value();
    let deposited_j_per_decay: f64 = assembly
        .model()
        .branches()
        .iter()
        .filter(|branch| branch.parent == parent)
        .map(|branch| branch.probability * branch.energy_per_decay.deposited_j)
        .sum();

    // Burnup accepts the original evaluated law and can add irradiation later.
    // This case explicitly supplies zero flux and no external isotope source.
    let groups = Arc::new(NeutronGroups::new(
        "unused-zero-flux-example-group".into(),
        vec![0.0, 1e-12],
    )?);
    let limits = BurnupTransportLimits::default();
    let model = Arc::new(BurnupModel::new(
        assembly.burnup_data(Arc::clone(&groups), 200.0, 1000.0, || Ok(()))?,
        limits.burnup,
    )?);
    let n = model.nuclides().len();
    let initial_nuclei = 1e12;
    let mut initial = vec![0.0; n];
    initial[parent] = initial_nuclei;
    let volumes_m3 = vec![1.0, 2.0];
    let volume_flow_m3_s = lambda * volumes_m3[0];
    let transport = NuclideTransport {
        transfers: vec![
            NuclideTransfer {
                from_compartment: 0,
                to_compartment: 1,
                rates_per_s: vec![volume_flow_m3_s / volumes_m3[0]; n],
            },
            NuclideTransfer {
                from_compartment: 1,
                to_compartment: 0,
                rates_per_s: vec![volume_flow_m3_s / volumes_m3[1]; n],
            },
        ],
        outflows: Vec::new(),
    };
    let mut state = BurnupTransportState::new(
        vec![
            BurnupInventory::new(Arc::clone(&model), initial)?,
            BurnupInventory::new(model, vec![0.0; n])?,
        ],
        volumes_m3,
        transport.clone(),
        limits,
        || Ok(()),
    )?;
    let flux = GroupFlux::new(groups, vec![0.0])?;
    let fluxes = [flux.clone(), flux];
    let supply = vec![vec![0.0; n]; 2];
    let budget = BurnupTransportBudget {
        limits,
        nuclei_tolerance: 0.1,
        energy_tolerance_j: 1e-12,
    };
    let step_s = half_life_s / 4.0;
    println!(
        "time_s,half_lives,h3_cell0,h3_cell1,he3_cell0,he3_cell1,deposited_cell0_j,deposited_cell1_j"
    );
    for sample in 0..=8 {
        if sample > 0 {
            (state, _) = state.trial_with_inputs(
                step_s,
                &fluxes,
                &[300.0, 300.0],
                &transport,
                &supply,
                budget,
                || Ok(()),
            )?;
        }
        let time_s = sample as f64 * step_s;
        let cells = state.inventories();
        // For common isotope transfer rates, spatial exchange and decay commute.
        // The unequal-volume equilibrium share is 1/3 in compartment zero.
        let equilibrium = 1.0 / 3.0;
        let exchange_rate = 1.5 * lambda;
        let fraction0 = equilibrium + (1.0 - equilibrium) * attenuate(1.0, exchange_rate * time_s);
        let parents = attenuate(initial_nuclei, lambda * time_s);
        let daughters = initial_nuclei * decayed_fraction(lambda * time_s);
        let decays0 = initial_nuclei
            * (equilibrium * decayed_fraction(lambda * time_s)
                + (1.0 - equilibrium) * lambda / (lambda + exchange_rate)
                    * decayed_fraction((lambda + exchange_rate) * time_s));
        for (cell, fraction, decays) in [
            (&cells[0], fraction0, decays0),
            (&cells[1], 1.0 - fraction0, daughters - decays0),
        ] {
            check(cell.populations()[parent], parents * fraction, 0.1)?;
            check(cell.populations()[daughter], daughters * fraction, 0.1)?;
            check(
                cell.cumulative_released().deposited_j,
                decays * deposited_j_per_decay,
                1e-14,
            )?;
        }
        println!(
            "{time_s:.9e},{:.2},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e}",
            sample as f64 / 4.0,
            cells[0].populations()[parent],
            cells[1].populations()[parent],
            cells[0].populations()[daughter],
            cells[1].populations()[daughter],
            cells[0].cumulative_released().deposited_j,
            cells[1].cumulative_released().deposited_j,
        );
    }
    Ok(())
}

fn check(actual: f64, expected: f64, allowance: f64) -> Result<(), Box<dyn std::error::Error>> {
    if !actual.is_finite() || (actual - expected).abs() > allowance {
        return Err(format!(
            "circulating-decay analytic comparison failed: actual={actual:e}, expected={expected:e}, allowance={allowance:e}"
        )
        .into());
    }
    Ok(())
}
