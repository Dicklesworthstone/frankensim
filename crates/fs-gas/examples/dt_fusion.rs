//! Run with: cargo run -p fs-gas --example dt_fusion
//! Constructed zero-dimensional ideal plasma; no confinement/ignition claim.

#[path = "dt_fusion/source.rs"]
mod source;

use fs_gas::fusion::{
    DtFusionHeating, DtFusionSolve, DtFusionTransport, advance_dt_fusion,
    advance_dt_fusion_with_transport,
};
use fs_qty::{Power, Temperature, Time};
use fs_thermochem::coulomb::{
    CoulombLimitsV1, CoulombModelV1, CoulombScreeningV1, CoulombShortRangeCutoffV1,
};
use fs_thermochem::fusion::KELVIN_PER_KEV;
use fs_thermochem::fusion::plasma::DtFusionPlasmaLimitsV1;
use fs_thermochem::plasma_radiation::{
    BremsstrahlungControlsV1, BremsstrahlungGauntFactorV1, ElectronIonBremsstrahlungV1,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let (model, initial) = source::plasma([5.0e19, 5.0e19, 0.0], 1.0, [10.0, 5.0]);
    let mut controls = DtFusionSolve::new(
        Temperature::new(0.2 * KELVIN_PER_KEV),
        Temperature::new(100.0 * KELVIN_PER_KEV),
        Temperature::new(1.0e4),
        Temperature::new(20.0 * KELVIN_PER_KEV),
        DtFusionPlasmaLimitsV1 {
            maximum_coulomb_coupling: 0.01,
            maximum_thermal_phase_space_density: 0.01,
            maximum_thermal_energy_rest_energy_ratio: 0.05,
        },
    );
    controls.relative_tolerance = 1.0e-8;
    controls.absolute_amount_tolerance_mol = 1.0e-20;
    controls.initial_step_s = 1.0e-3;
    controls.maximum_step_s = 0.1;
    let heating = DtFusionHeating {
        alpha_electron_fraction: 0.8,
        electron_power: Power::new(0.0),
        heavy_power: Power::new(0.0),
    };
    let isolated = advance_dt_fusion(&initial, &model, heating, Time::new(0.1), controls, || {
        false
    })?;
    let transport = DtFusionTransport {
        coulomb: Some(CoulombModelV1::new_with_cutoff(
            &initial,
            CoulombScreeningV1::StaticAllChargedDebye,
            CoulombShortRangeCutoffV1::QuantumDiffraction,
            CoulombLimitsV1 {
                minimum_coulomb_logarithm: 4.0,
                minimum_particles_per_debye_sphere: 10.0,
                maximum_coulomb_coupling: 0.01,
                maximum_thermal_phase_space_density: 0.01,
                maximum_thermal_energy_rest_energy_ratio: 0.05,
                maximum_quantum_to_classical_ratio: 0.1,
                maximum_relative_drift_thermal_ratio: 0.1,
            },
            || false,
        )?),
        bremsstrahlung: Some(ElectronIonBremsstrahlungV1::new(
            model.thermodynamics().clone(),
            BremsstrahlungControlsV1 {
                frequency_integrated_gaunt_factor: BremsstrahlungGauntFactorV1::ItohNonrelativistic,
                minimum_electron_temperature_k: 1e6,
                maximum_electron_temperature_k: 20.0 * KELVIN_PER_KEV,
                maximum_electron_thermal_rest_energy_ratio: 0.05,
                maximum_ion_to_electron_thermal_variance_ratio: 0.05,
                maximum_electron_thermal_phase_space_density: 0.01,
            },
        )?),
    };
    let result = advance_dt_fusion_with_transport(
        &initial,
        &model,
        heating,
        &transport,
        Time::new(0.1),
        controls,
        || false,
    )?;
    let state = result.state();
    let balance = result.balance();
    println!("D-T burn: rigid 1 m3 parcel, nD = nT = 5e19 /m3, duration 0.1 s");
    println!(
        "Closure: instantaneous local alpha thermalization (80% electrons); all neutrons and photons escape."
    );
    println!("Source: ideal monatomic calorics and bounded Maxwellian Bosch-Hale rates.");
    println!(
        "Exchange: static all-carrier Debye screening with an explicit quantum-diffraction cutoff."
    );
    println!(
        "Radiation: electron-ion free-free emission with the bounded Itoh temperature/charge Gaunt fit."
    );
    println!(
        "No confinement, fast-particle slowing, opacity or other emission channels are modeled."
    );
    println!(
        "Without exchange/photon losses: ion/electron kT {:.9} / {:.9} keV; extent {:.12e} mol",
        isolated.state().heavy_temperature().value() / KELVIN_PER_KEV,
        isolated.state().electron_temperature().value() / KELVIN_PER_KEV,
        isolated.balance().reaction_extent_mol
    );
    println!(
        "Final ion/electron kT: {:.9} / {:.9} keV",
        state.heavy_temperature().value() / KELVIN_PER_KEV,
        state.electron_temperature().value() / KELVIN_PER_KEV
    );
    println!(
        "D and T consumed; He4 and neutrons created: {:.12e} mol each",
        balance.reaction_extent_mol
    );
    println!(
        "Rest-mass energy released: {:.9} J",
        balance.nuclear_energy_released_j
    );
    println!(
        "Alpha energy deposited: {:.9} J",
        balance.alpha_energy_deposited_j
    );
    println!(
        "Reacting-ion energy removed: {:.9} J",
        balance.reacting_ion_energy_removed_j
    );
    println!(
        "Neutron energy escaped: {:.9} J",
        balance.neutron_energy_escaped_j
    );
    println!(
        "Signed electron-to-ion Coulomb transfer: {:.9} J",
        balance.electron_to_heavy_coulomb_energy_j
    );
    println!(
        "Photon energy escaped: {:.9} J",
        balance.photon_energy_escaped_j
    );
    println!("Final pressure: {:.9} Pa", state.pressure().value());
    println!(
        "Energy-account residual: {:.6e} J",
        balance.total_energy_residual_j
    );
    println!(
        "Accepted/rejected intervals: {} / {}",
        balance.accepted_steps, balance.rejected_steps
    );
    Ok(())
}
