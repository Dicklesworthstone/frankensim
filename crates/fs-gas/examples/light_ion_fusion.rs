//! Primary and secondary light-ion fusion on one original material state.
//! Run with: cargo run -p fs-gas --example light_ion_fusion

#[path = "light_ion_fusion/source.rs"]
mod source;

use fs_gas::fusion::network::{
    FusionNetworkChannels, FusionNetworkHeating, FusionNetworkResult, FusionNetworkSolve,
    FusionNetworkTransport, advance_fusion_network,
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

fn controls() -> FusionNetworkSolve {
    let mut c = FusionNetworkSolve::new(
        Temperature::new(0.5 * KELVIN_PER_KEV),
        Temperature::new(100.0 * KELVIN_PER_KEV),
        Temperature::new(1e6),
        Temperature::new(25.0 * KELVIN_PER_KEV),
        DtFusionPlasmaLimitsV1 {
            maximum_coulomb_coupling: 0.01,
            maximum_thermal_phase_space_density: 0.01,
            maximum_thermal_energy_rest_energy_ratio: 0.05,
        },
    );
    c.relative_tolerance = 1e-8;
    c.absolute_amount_tolerance_mol = 1e-20;
    c.initial_step_s = 1e-3;
    c.maximum_step_s = 0.1;
    c
}

fn report(label: &str, result: &FusionNetworkResult) {
    let b = result.balance();
    let state = result.state();
    println!("{label}");
    println!(
        "  final ion/electron kT: {:.9} / {:.9} keV",
        state.heavy_temperature().value() / KELVIN_PER_KEV,
        state.electron_temperature().value() / KELVIN_PER_KEV
    );
    for (label, index) in ["p", "D", "T", "He3", "He4"]
        .into_iter()
        .zip(result.plasma().indices())
    {
        let amount = state.amounts()[index].value();
        println!("  final {label}: {amount:.12e} mol");
    }
    for (reaction, extent) in ["DT", "DD->He3+n", "DD->T+p", "DHe3"]
        .into_iter()
        .zip(b.reaction_extents_mol)
    {
        println!("  {reaction} events: {extent:.12e} mol");
    }
    println!("  nuclear Q: {:.9} J", b.nuclear_energy_released_j);
    println!(
        "  charged-product electron/ion deposition: {:.9} / {:.9} J",
        b.charged_product_electron_energy_j, b.charged_product_heavy_energy_j
    );
    println!(
        "  reacting-ion thermal debit: {:.9} J",
        b.reacting_ion_energy_removed_j
    );
    println!(
        "  escaped neutron amount/energy: {:.12e} mol / {:.9} J",
        b.neutron_amount_mol, b.neutron_energy_escaped_j
    );
    println!(
        "  signed electron-to-ion Coulomb transfer: {:.9} J",
        b.electron_to_heavy_coulomb_energy_j
    );
    println!("  escaped photons: {:.9} J", b.photon_energy_escaped_j);
    println!(
        "  energy residual: {:.6e} J; baryon/charge residual: {:.3e} / {:.3e}",
        b.total_energy_residual_j, b.relative_baryon_residual, b.relative_charge_residual
    );
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("Light-ion fusion: original ideal-monatomic sources and bounded Maxwellian rates.");
    println!(
        "Closure: homogeneous rigid parcel, immediate charged-product thermalization, all neutrons escape."
    );
    println!(
        "Selected radiation assumes all emitted photons escape; no confinement or fast-particle slowing."
    );

    let (plasma, initial) = source::plasma([0.0, 2e21, 0.0, 0.0, 0.0], 1.0, [20.0, 8.0]);
    let heating =
        FusionNetworkHeating::uniform_charged_fraction(0.8, Power::new(0.0), Power::new(0.0));
    let primary = advance_fusion_network(
        &initial,
        &plasma,
        FusionNetworkChannels::new([false, true, true, false]),
        heating,
        &FusionNetworkTransport::default(),
        Time::new(1.0),
        controls(),
        || false,
    )?;
    let complete = advance_fusion_network(
        &initial,
        &plasma,
        FusionNetworkChannels::ALL,
        heating,
        &FusionNetworkTransport::default(),
        Time::new(1.0),
        controls(),
        || false,
    )?;
    println!(
        "Case A: pure D, nD=2e21/m3, V=1m3, Ti/Te=20/8keV, 1s, 80% charged deposition to electrons."
    );
    println!("Thermal operators are explicitly absent in both A comparisons.");
    report("Primary DD channels only:", &primary);
    report("All four channels, with no T or He3 seed:", &complete);

    let (plasma, initial) = source::plasma([2e19, 3e20, 7e19, 4e19, 1e19], 0.025, [15.0, 4.0]);
    let transport = FusionNetworkTransport {
        particles: None,
        energy_confinement: None,
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
            plasma.thermodynamics().clone(),
            BremsstrahlungControlsV1 {
                frequency_integrated_gaunt_factor: BremsstrahlungGauntFactorV1::ItohNonrelativistic,
                minimum_electron_temperature_k: 1e6,
                maximum_electron_temperature_k: 25.0 * KELVIN_PER_KEV,
                maximum_electron_thermal_rest_energy_ratio: 0.05,
                maximum_ion_to_electron_thermal_variance_ratio: 0.05,
                maximum_electron_thermal_phase_space_density: 0.01,
            },
        )?),
    };
    let heating = FusionNetworkHeating {
        product_electron_fractions: [[0.8, 0.0], [0.65, 0.0], [0.4, 0.9], [0.75, 0.9]],
        electron_power: Power::new(0.0),
        heavy_power: Power::new(0.0),
    };
    let coupled = advance_fusion_network(
        &initial,
        &plasma,
        FusionNetworkChannels::ALL,
        heating,
        &transport,
        Time::new(0.2),
        controls(),
        || false,
    )?;
    println!(
        "Case B: n[p,D,T,He3,He4]=[2e19,3e20,7e19,4e19,1e19]/m3, V=.025m3, Ti/Te=15/4keV, .2s."
    );
    println!(
        "Quantum-cutoff Coulomb exchange and bounded Itoh Gaunt radiation; authored per-product deposition fractions."
    );
    report("All four channels with thermal feedback:", &coupled);
    Ok(())
}
