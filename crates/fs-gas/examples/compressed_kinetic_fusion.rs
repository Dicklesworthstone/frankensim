//! Finite-time alpha heating during prescribed compression, hold and expansion.
//! Run with: cargo run -p fs-gas --example compressed_kinetic_fusion
//!
//! Both cases retain the original ideal-monatomic material sources and the
//! same Coulomb/radiation models. The moving case assumes a confined isotropic
//! parcel; the volume history supplies boundary motion rather than solving
//! piston or confinement dynamics. Alpha cell amounts move without remapping
//! under geometry alone, while their physical energies obey E*V^(2/3)=constant.
//! Energy-grid accuracy remains separate from adaptive time-integration error.

#[path = "light_ion_fusion/source.rs"]
mod source;

use fs_gas::fusion::kinetic::{
    AlphaEnergyGrid, KineticAlphaTransport, KineticFusionResult, KineticFusionState,
    advance_kinetic_fusion, advance_kinetic_fusion_with_volume,
};
use fs_gas::fusion::network::{
    FusionNetworkChannels, FusionNetworkHeating, FusionNetworkSolve, FusionNetworkTransport,
};
use fs_gas::fusion::volume::{FusionVolumeHistory, FusionVolumeKnot};
use fs_qty::{Energy, Power, Temperature, Time, Volume};
use fs_thermochem::coulomb::test_particle::TestParticleCoulombV1;
use fs_thermochem::coulomb::{
    CoulombLimitsV1, CoulombModelV1, CoulombScreeningV1, CoulombShortRangeCutoffV1,
};
use fs_thermochem::electron_collision::AVOGADRO_PER_MOL;
use fs_thermochem::fusion::KELVIN_PER_KEV;
use fs_thermochem::fusion::plasma::DtFusionPlasmaLimitsV1;
use fs_thermochem::plasma_radiation::{
    BremsstrahlungControlsV1, BremsstrahlungGauntFactorV1, ElectronIonBremsstrahlungV1,
};

const KEV_J: f64 = 1.602_176_634e-16;

fn controls() -> FusionNetworkSolve {
    let mut controls = FusionNetworkSolve::new(
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
    controls.relative_tolerance = 1e-5;
    controls.absolute_amount_tolerance_mol = 1e-18;
    controls.initial_step_s = 1e-4;
    controls.maximum_step_s = 5e-3;
    controls
}

fn report(label: &str, result: &KineticFusionResult) {
    let state = result.state();
    let balance = result.balance();
    let alpha = result.alpha_balance();
    println!("{label}");
    println!(
        "  elapsed time: {:.6} s; final volume: {:.9} m3",
        result.elapsed_time().value(),
        state.volume().value()
    );
    match result.volume_history() {
        Some(history) => {
            println!("  retained log-volume history [time s, volume m3]:");
            for knot in history.knots() {
                println!("    {:.6}, {:.9}", knot.time.value(), knot.volume.value());
            }
        }
        None => println!("  retained geometry: rigid, constant input volume"),
    }
    println!(
        "  final background ion/electron kT: {:.9} / {:.9} keV",
        state.background_heavy_temperature().value() / KELVIN_PER_KEV,
        state.electron_temperature().value() / KELVIN_PER_KEV
    );
    println!(
        "  thermal ion/electron stores: {:.9} / {:.9} J; kinetic alpha store: {:.9} J",
        state.thermal_heavy_energy().value(),
        state.thermal_electron_energy().value(),
        state.alpha_energy().value()
    );
    println!(
        "  total stored-energy change: {:.9} J",
        state.internal_energy().value() - result.initial_state().internal_energy().value()
    );
    let pressure_factor = (2.0 / 3.0) / state.volume().value();
    println!(
        "  ion/electron/alpha pressures: {:.9} / {:.9} / {:.9} Pa; total: {:.9} Pa",
        pressure_factor * state.thermal_heavy_energy().value(),
        pressure_factor * state.thermal_electron_energy().value(),
        pressure_factor * state.alpha_energy().value(),
        state.pressure().value()
    );
    println!(
        "  retained kinetic alphas: {:.12e} particles ({:.12e} mol)",
        state.alpha_amount().value() * AVOGADRO_PER_MOL,
        state.alpha_amount().value()
    );
    let physical_grid = state.energy_grid();
    println!(
        "  published physical alpha grid: {} centers, {:.9} to {:.9} keV",
        physical_grid.len(),
        physical_grid.energies()[0].value() / KEV_J,
        physical_grid.energies()[physical_grid.len() - 1].value() / KEV_J
    );
    for (index, channel) in ["DT", "DD->He3+n", "DD->T+p", "DHe3"]
        .into_iter()
        .enumerate()
    {
        let products = balance.product_kinetic_energies_j[index];
        println!(
            "  {channel}: {:.12e} mol events; product kinetic energies {:.9} / {:.9} J",
            balance.reaction_extents_mol[index], products[0], products[1]
        );
    }
    println!("  nuclear Q: {:.9} J", balance.nuclear_energy_released_j);
    println!(
        "  reacting-ion thermal energy removed: {:.9} J",
        balance.reacting_ion_energy_removed_j
    );
    println!(
        "  alpha births: {:.12e} mol / {:.9} J",
        alpha.born_alpha_amount_mol, alpha.born_alpha_energy_j
    );
    println!(
        "  signed alpha collision heat to ions/electrons: {:.9} / {:.9} J",
        alpha.deposited_heavy_energy_j, alpha.deposited_electron_energy_j
    );
    println!(
        "  prompt non-alpha heat to ions/electrons: {:.9} / {:.9} J",
        balance.charged_product_heavy_energy_j, balance.charged_product_electron_energy_j
    );
    println!(
        "  signed electron-to-ion background Coulomb transfer: {:.9} J",
        balance.electron_to_heavy_coulomb_energy_j
    );
    println!(
        "  escaped neutron amount/energy: {:.12e} mol / {:.9} J; photons: {:.9} J",
        balance.neutron_amount_mol,
        balance.neutron_energy_escaped_j,
        balance.photon_energy_escaped_j
    );
    println!(
        "  external heat to ions/electrons: {:.9} / {:.9} J",
        balance.external_heavy_energy_j, balance.external_electron_energy_j
    );
    println!(
        "  signed boundary work on ions/electrons/alphas: {:.9} / {:.9} / {:.9} J",
        balance.boundary_work_on_heavy_j,
        balance.boundary_work_on_electrons_j,
        alpha.boundary_work_on_alphas_j
    );
    println!(
        "  total mechanical work on the parcel: {:.9} J",
        result.total_boundary_work().value()
    );
    println!(
        "  alpha residual (delta Ualpha + collision heat - births - work): {:.6e} J",
        alpha.alpha_energy_residual_j
    );
    println!(
        "  original-source ion/electron energy residuals: {:.6e} / {:.6e} J",
        balance.heavy_energy_residual_j, balance.electron_energy_residual_j
    );
    println!(
        "  full energy residual (delta U + neutrons + photons - Q - heat - all work): {:.6e} J",
        balance.total_energy_residual_j
    );
    println!(
        "  relative baryon/charge residuals: {:.3e} / {:.3e}; species/alpha residuals: {:.3e} / {:.3e} mol",
        balance.relative_baryon_residual,
        balance.relative_charge_residual,
        balance.maximum_stoichiometry_residual_mol,
        alpha.alpha_amount_residual_mol
    );
    println!(
        "  uppermost alpha-bin fraction: {:.6e}",
        alpha.upper_bin_fraction
    );
    println!(
        "  accepted/rejected steps: {} / {}; largest normalized time-error estimate: {:.6}",
        balance.accepted_steps, balance.rejected_steps, balance.maximum_normalized_local_error
    );
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let (plasma, initial) = source::plasma([0.0, 5e19, 5e19, 0.0, 0.0], 1.0, [10.0, 5.0]);
    let controls = controls();
    let channels = FusionNetworkChannels::new([true, false, false, false]);
    let heating =
        FusionNetworkHeating::uniform_charged_fraction(0.0, Power::new(0.0), Power::new(0.0));
    let grid = AlphaEnergyGrid::uniform_speed(
        Energy::new(0.05 * KEV_J),
        Energy::new(5000.0 * KEV_J),
        256,
        || false,
    )?;
    let kinetic_initial =
        KineticFusionState::from_thermal(&initial, &plasma, &grid, controls.plasma_limits, || {
            false
        })?;
    let background = CoulombModelV1::new_with_cutoff(
        &initial,
        CoulombScreeningV1::ElectronDebye,
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
    )?;
    let alpha_transport = KineticAlphaTransport::new(
        TestParticleCoulombV1::new(background.clone(), plasma.indices()[4], 0.01)?,
        0.05,
        0.25,
        1e-4,
        100_000,
    )?;
    let thermal_transport = FusionNetworkTransport {
        coulomb: Some(background),
        bremsstrahlung: Some(ElectronIonBremsstrahlungV1::new(
            plasma.thermodynamics().clone(),
            BremsstrahlungControlsV1 {
                frequency_integrated_gaunt_factor: BremsstrahlungGauntFactorV1::Constant(1.2),
                minimum_electron_temperature_k: controls.minimum_electron_temperature_k,
                maximum_electron_temperature_k: controls.maximum_electron_temperature_k,
                maximum_electron_thermal_rest_energy_ratio: 0.05,
                maximum_ion_to_electron_thermal_variance_ratio: 0.1,
                maximum_electron_thermal_phase_space_density: 0.01,
            },
        )?),
    };
    let history = FusionVolumeHistory::new(
        [(0.0, 1.0), (0.04, 0.6), (0.06, 0.6), (0.1, 1.0)]
            .into_iter()
            .map(|(time, volume)| FusionVolumeKnot {
                time: Time::new(time),
                volume: Volume::new(volume),
            })
            .collect(),
        || false,
    )?;

    println!("D–T only: nD=nT=5e19/m3, V=1m3, initial Ti/Te=10/5keV, duration=0.1s.");
    println!(
        "Both cases start with zero He4 and use the same original material and collision sources."
    );
    println!(
        "Kinetic alpha drift/diffusion, background Coulomb exchange, Gaunt1.2 free-free radiation."
    );
    println!(
        "Electron-Debye screening and thermal-reference quantum-diffraction cutoff in both cases."
    );
    println!("No prompt alpha heat or external heat; neutrons and optically thin photons escape.");
    println!(
        "The moving case assumes a confined isotropic parcel with prescribed log-volume motion."
    );
    println!("No piston dynamics, confinement solution or anisotropic magnetic work is inferred.");
    println!(
        "Initial alpha grid: 256 centers uniform in speed, 0.05–5000keV, reflecting energy faces."
    );
    println!(
        "Geometry scales physical centers and widths with V^(-2/3), without remapping bin amounts."
    );
    println!(
        "Collisions and births use physical energies; published grids contain those physical energies."
    );
    println!(
        "Grid accuracy and mean-energy birth spectra remain separate from the adaptive time tolerance."
    );
    println!(
        "Alpha/thermal-ion number <= .05, alpha/thermal energy <= .25, upper-bin fraction <= 1e-4."
    );
    println!(
        "Mechanical work is positive on compression and negative on expansion in each component."
    );

    let rigid = advance_kinetic_fusion(
        &kinetic_initial,
        channels,
        heating,
        &thermal_transport,
        &alpha_transport,
        history.duration(),
        controls,
        || false,
    )?;
    let moving = advance_kinetic_fusion_with_volume(
        &kinetic_initial,
        channels,
        heating,
        &thermal_transport,
        &alpha_transport,
        &history,
        controls,
        || false,
    )?;
    report("Rigid kinetic baseline:", &rigid);
    report("Kinetic compression, hold and expansion:", &moving);
    Ok(())
}
