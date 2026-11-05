//! Coupled fusion during prescribed compression, a hold, and expansion.
//! Run with: cargo run -p fs-gas --example compressed_fusion
//!
//! Illustrative homogeneous Maxwellian parcel with original ideal-monatomic
//! source cards. Geometry is authored, not inferred from confinement or piston
//! dynamics. Charged products thermalize immediately; neutrons and photons
//! escape. This is not a fast-particle or radiation-transport calculation.

#[path = "light_ion_fusion/source.rs"]
mod source;

use fs_gas::fusion::network::{
    FusionNetworkChannels, FusionNetworkHeating, FusionNetworkResult, FusionNetworkSolve,
    FusionNetworkTransport, advance_fusion_network, advance_fusion_network_with_volume,
};
use fs_gas::fusion::volume::{FusionVolumeHistory, FusionVolumeKnot};
use fs_qty::{Power, Temperature, Time, Volume};
use fs_thermochem::coulomb::{
    CoulombLimitsV1, CoulombModelV1, CoulombScreeningV1, CoulombShortRangeCutoffV1,
};
use fs_thermochem::fusion::KELVIN_PER_KEV;
use fs_thermochem::fusion::plasma::DtFusionPlasmaLimitsV1;
use fs_thermochem::plasma_radiation::{
    BremsstrahlungControlsV1, BremsstrahlungGauntFactorV1, ElectronIonBremsstrahlungV1,
};

fn controls() -> FusionNetworkSolve {
    let mut controls = FusionNetworkSolve::new(
        Temperature::new(0.5 * KELVIN_PER_KEV),
        Temperature::new(100.0 * KELVIN_PER_KEV),
        Temperature::new(1e6),
        Temperature::new(20.0 * KELVIN_PER_KEV),
        DtFusionPlasmaLimitsV1 {
            maximum_coulomb_coupling: 0.01,
            maximum_thermal_phase_space_density: 0.01,
            maximum_thermal_energy_rest_energy_ratio: 0.05,
        },
    );
    controls.relative_tolerance = 1e-8;
    controls.absolute_amount_tolerance_mol = 1e-20;
    controls.initial_step_s = 1e-3;
    controls.maximum_step_s = 0.1;
    controls
}

fn report(label: &str, result: &FusionNetworkResult) {
    let balance = result.balance();
    let state = result.state();
    println!("{label}");
    println!(
        "  elapsed time: {:.6} s; final volume: {:.9} m3",
        result.elapsed_time().value(),
        state.volume().value()
    );
    println!(
        "  final ion/electron kT: {:.9} / {:.9} keV",
        state.heavy_temperature().value() / KELVIN_PER_KEV,
        state.electron_temperature().value() / KELVIN_PER_KEV
    );
    for (channel, extent) in ["DT", "DD->He3+n", "DD->T+p", "DHe3"]
        .into_iter()
        .zip(balance.reaction_extents_mol)
    {
        println!("  {channel} events: {extent:.12e} mol");
    }
    println!("  nuclear Q: {:.9} J", balance.nuclear_energy_released_j);
    println!(
        "  escaped neutron amount/energy: {:.12e} mol / {:.9} J",
        balance.neutron_amount_mol, balance.neutron_energy_escaped_j
    );
    println!(
        "  escaped photons: {:.9} J",
        balance.photon_energy_escaped_j
    );
    println!(
        "  signed electron-to-ion Coulomb transfer: {:.9} J",
        balance.electron_to_heavy_coulomb_energy_j
    );
    println!(
        "  boundary work on ions/electrons: {:.9} / {:.9} J",
        balance.boundary_work_on_heavy_j, balance.boundary_work_on_electrons_j
    );
    println!(
        "  energy residual (delta U + neutrons + photons - Q - heat - work): {:.6e} J",
        balance.total_energy_residual_j
    );
    println!(
        "  relative baryon/charge residual: {:.3e} / {:.3e}; species residual: {:.3e} mol",
        balance.relative_baryon_residual,
        balance.relative_charge_residual,
        balance.maximum_stoichiometry_residual_mol
    );
    println!(
        "  accepted/rejected steps: {} / {}; largest normalized local error: {:.6}",
        balance.accepted_steps, balance.rejected_steps, balance.maximum_normalized_local_error
    );
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let (plasma, initial) = source::plasma([2e19, 3e20, 7e19, 4e19, 1e19], 0.025, [15.0, 4.0]);
    let heating = FusionNetworkHeating {
        product_electron_fractions: [[0.8, 0.0], [0.65, 0.0], [0.4, 0.9], [0.75, 0.9]],
        electron_power: Power::new(0.0),
        heavy_power: Power::new(0.0),
    };
    let transport = FusionNetworkTransport {
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
                frequency_integrated_gaunt_factor: BremsstrahlungGauntFactorV1::Constant(1.2),
                minimum_electron_temperature_k: 1e6,
                maximum_electron_temperature_k: 20.0 * KELVIN_PER_KEV,
                maximum_electron_thermal_rest_energy_ratio: 0.05,
                maximum_ion_to_electron_thermal_variance_ratio: 0.05,
                maximum_electron_thermal_phase_space_density: 0.01,
            },
        )?),
    };
    let history = FusionVolumeHistory::new(
        [(0.0, 0.025), (0.04, 0.0125), (0.06, 0.0125), (0.1, 0.02)]
            .into_iter()
            .map(|(time, volume)| FusionVolumeKnot {
                time: Time::new(time),
                volume: Volume::new(volume),
            })
            .collect(),
        || false,
    )?;

    println!("Homogeneous fusion with prescribed log-linear volume and two thermal pools.");
    println!("Initial n[p,D,T,He3,He4]=[2e19,3e20,7e19,4e19,1e19]/m3, V=.025m3, Ti/Te=15/4keV.");
    println!(
        "All four channels, quantum-cutoff Coulomb exchange, constant Gaunt 1.2, zero external heat."
    );
    println!(
        "Electron deposition fractions by product: DT [.8,0], DDn [.65,0], DDp [.4,.9], DHe3 [.75,.9]."
    );
    println!(
        "Illustrative ideal-monatomic sources; immediate charged-product thermalization; all neutrons and photons escape."
    );
    println!("No inferred confinement, piston dynamics, or fast-particle slowing.");
    println!("Prescribed time [s] / volume [m3] knots:");
    for knot in history.knots() {
        println!("  {:.3} / {:.6}", knot.time.value(), knot.volume.value());
    }
    println!(
        "Work is positive onto each plasma pool during compression and negative during expansion."
    );

    let rigid = advance_fusion_network(
        &initial,
        &plasma,
        FusionNetworkChannels::ALL,
        heating,
        &transport,
        history.duration(),
        controls(),
        || false,
    )?;
    let compressed = advance_fusion_network_with_volume(
        &initial,
        &plasma,
        FusionNetworkChannels::ALL,
        heating,
        &transport,
        &history,
        controls(),
        || false,
    )?;
    report("Same-duration rigid baseline:", &rigid);
    report("Prescribed compression, hold and expansion:", &compressed);
    Ok(())
}
