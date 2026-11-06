//! D–T burn with finite-time, isotropic alpha slowing and energy diffusion.
//! Run with: cargo run -p fs-gas --example kinetic_fusion
//!
//! Original ideal-monatomic sources, a rigid homogeneous parcel and escaping
//! neutrons. Only He4 is kinetic; fuels and electrons remain Maxwellian. The
//! finite reflecting energy grid is a numerical approximation, and grid
//! convergence must be assessed separately from adaptive time-step convergence.

#[path = "light_ion_fusion/source.rs"]
mod source;

use fs_gas::fusion::kinetic::{
    AlphaEnergyGrid, KineticAlphaTransport, KineticFusionResult, KineticFusionState,
    advance_kinetic_fusion,
};
use fs_gas::fusion::network::{
    FusionNetworkBalance, FusionNetworkChannels, FusionNetworkHeating, FusionNetworkResult,
    FusionNetworkSolve, FusionNetworkTransport, advance_fusion_network,
};
use fs_qty::{Energy, Power, Temperature, Time};
use fs_thermochem::coulomb::test_particle::TestParticleCoulombV1;
use fs_thermochem::coulomb::{
    CoulombLimitsV1, CoulombModelV1, CoulombScreeningV1, CoulombShortRangeCutoffV1,
};
use fs_thermochem::electron_collision::{AVOGADRO_PER_MOL, BOLTZMANN_J_PER_K};
use fs_thermochem::fusion::KELVIN_PER_KEV;
use fs_thermochem::fusion::plasma::DtFusionPlasmaLimitsV1;

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

fn report_balance(balance: FusionNetworkBalance) {
    println!("  DT events: {:.12e} mol", balance.reaction_extents_mol[0]);
    println!("  nuclear Q: {:.9} J", balance.nuclear_energy_released_j);
    println!(
        "  escaped neutron amount/energy: {:.12e} mol / {:.9} J",
        balance.neutron_amount_mol, balance.neutron_energy_escaped_j
    );
    println!(
        "  reacting-ion thermal energy removed: {:.9} J",
        balance.reacting_ion_energy_removed_j
    );
    println!(
        "  total energy residual (delta thermal + delta kinetic + neutrons - Q): {:.6e} J",
        balance.total_energy_residual_j
    );
    println!(
        "  ion/electron energy residuals: {:.6e} / {:.6e} J",
        balance.heavy_energy_residual_j, balance.electron_energy_residual_j
    );
    println!(
        "  relative baryon/charge residuals: {:.3e} / {:.3e}; species residual: {:.3e} mol",
        balance.relative_baryon_residual,
        balance.relative_charge_residual,
        balance.maximum_stoichiometry_residual_mol
    );
    println!(
        "  accepted/rejected steps: {} / {}; largest normalized local error: {:.6}",
        balance.accepted_steps, balance.rejected_steps, balance.maximum_normalized_local_error
    );
}

fn report_prompt(result: &FusionNetworkResult) {
    let state = result.state();
    let balance = result.balance();
    let alpha_amount = state.amounts()[result.plasma().indices()[4]].value();
    let alpha_number = alpha_amount * AVOGADRO_PER_MOL;
    let alpha_pressure = (alpha_number / state.volume().value())
        * BOLTZMANN_J_PER_K
        * state.heavy_temperature().value();
    println!("Prompt baseline: explicit 80% alpha birth energy to electrons, 20% to ions.");
    println!(
        "  final ion/electron kT: {:.9} / {:.9} keV",
        state.heavy_temperature().value() / KELVIN_PER_KEV,
        state.electron_temperature().value() / KELVIN_PER_KEV
    );
    println!(
        "  thermal ion/electron stores: {:.9} / {:.9} J; separate kinetic alpha store: 0 J",
        state.heavy_internal_energy().value(),
        state.electron_internal_energy().value()
    );
    println!("  retained thermal alphas: {alpha_number:.12e} particles ({alpha_amount:.12e} mol)");
    println!(
        "  alpha/total pressure: {alpha_pressure:.9} / {:.9} Pa",
        state.pressure().value()
    );
    println!(
        "  alpha heat to ions/electrons: {:.9} / {:.9} J",
        balance.charged_product_heavy_energy_j, balance.charged_product_electron_energy_j
    );
    report_balance(balance);
}

fn report_kinetic(result: &KineticFusionResult) {
    let state = result.state();
    let alpha = result.alpha_balance();
    let alpha_number = state.alpha_amount().value() * AVOGADRO_PER_MOL;
    let alpha_pressure = (2.0 / 3.0) * state.alpha_energy().value() / state.volume().value();
    println!(
        "Kinetic alpha evolution: collision-derived partition, with zero prompt alpha deposition."
    );
    println!(
        "  final background ion/electron kT: {:.9} / {:.9} keV",
        state.background_heavy_temperature().value() / KELVIN_PER_KEV,
        state.electron_temperature().value() / KELVIN_PER_KEV
    );
    println!(
        "  thermal ion/electron stores: {:.9} / {:.9} J; retained kinetic alpha store: {:.9} J",
        state.thermal_heavy_energy().value(),
        state.thermal_electron_energy().value(),
        state.alpha_energy().value()
    );
    println!(
        "  retained kinetic alphas: {alpha_number:.12e} particles ({:.12e} mol)",
        state.alpha_amount().value()
    );
    println!(
        "  alpha/total pressure: {alpha_pressure:.9} / {:.9} Pa",
        state.pressure().value()
    );
    println!(
        "  born alpha kinetic energy: {:.9} J; alpha heat to ions/electrons: {:.9} / {:.9} J",
        alpha.born_alpha_energy_j,
        alpha.deposited_heavy_energy_j,
        alpha.deposited_electron_energy_j
    );
    println!(
        "  alpha storage+heat-birth residual: {:.6e} J; alpha amount residual: {:.6e} mol",
        alpha.alpha_energy_residual_j, alpha.alpha_amount_residual_mol
    );
    println!("  uppermost-bin fraction: {:.6e}", alpha.upper_bin_fraction);
    report_balance(result.balance());
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let (plasma, initial) = source::plasma([0.0, 5e19, 5e19, 0.0, 0.0], 1.0, [10.0, 5.0]);
    let controls = controls();
    let channels = FusionNetworkChannels::new([true, false, false, false]);
    let duration = Time::new(0.1);
    let thermal_transport = FusionNetworkTransport::default();
    let prompt_heating =
        FusionNetworkHeating::uniform_charged_fraction(0.8, Power::new(0.0), Power::new(0.0));
    let kinetic_heating =
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
    let collisions = TestParticleCoulombV1::new(
        CoulombModelV1::new_with_cutoff(
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
        )?,
        plasma.indices()[4],
        0.01,
    )?;
    let alpha_transport = KineticAlphaTransport::new(collisions, 0.01, 0.20, 1e-6, 100_000)?;

    println!("D–T only: nD=nT=5e19/m3, V=1m3, initial Ti/Te=10/5keV, duration=0.1s.");
    println!(
        "Both cases start with exactly zero alpha inventory and share the original material sources."
    );
    println!(
        "Rigid homogeneous parcel; neutron escape; zero external heat and no other thermal transport."
    );
    println!(
        "Kinetic case: electron-Debye screening, thermal-reference quantum-diffraction cutoff."
    );
    println!("Alpha energy drift and diffusion remain active through the thermal tail.");
    println!("256 centers uniform in speed, 0.05–5000keV; reflecting energy boundaries.");
    println!(
        "Finite-grid accuracy needs separate refinement; the time-step tolerance is not a grid-error bound."
    );
    println!("Dilute limits: alpha/thermal-ion number <= .01, alpha/thermal energy <= .20.");
    println!(
        "Only alpha energy is kinetic; no pitch-angle, confinement, fast-secondary reaction or neutron transport."
    );

    let prompt = advance_fusion_network(
        &initial,
        &plasma,
        channels,
        prompt_heating,
        &thermal_transport,
        duration,
        controls,
        || false,
    )?;
    report_prompt(&prompt);
    let kinetic = advance_kinetic_fusion(
        &kinetic_initial,
        channels,
        kinetic_heating,
        &thermal_transport,
        &alpha_transport,
        duration,
        controls,
        || false,
    )?;
    report_kinetic(&kinetic);
    Ok(())
}
