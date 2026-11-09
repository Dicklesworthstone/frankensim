//! Heat-driven vapor mass, radial dynamics and causal water-borne sound.
//!
//! Run: cargo run -p fs-couple --features cavitation --example vapor_cavitation
//! Compare zero heat with an explicitly authored 1 mW smooth heat pulse under
//! the same pressure cycle. The pulse is an external thermal port, not an
//! inferred liquid-conduction or finite-rate evaporation correlation.
//! Water source, phase density, latent calorics and acoustic propagation use
//! the original IF97 material. The native approximation is locally saturated,
//! spatially uniform pure vapor in an incompressible frozen bulk liquid.

use fs_aeroac::monopole::MonopoleLimits;
use fs_couple::cavitation::bubble::vapor::MaterialVaporBubble;
use fs_flux::cavitation::vapor::{
    VaporBubbleError, VaporBubbleLimits, VaporBubbleNumerics, VaporBubbleState,
    VaporBubbleThermalState,
};
use fs_flux::cavitation::{BubbleLimits, BubbleNumerics};
use fs_material::cavitation::SaturatedWaterCavitation;
use fs_qty::Temperature;
use fs_thermochem::water::if97::If97SaturationPair;

const RADIUS: f64 = 100e-6;
const TEMPERATURE: f64 = 373.15;
const DURATION: f64 = 40e-6;
const SAMPLES: usize = 160;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let pair = If97SaturationPair::at_temperature(Temperature::new(TEMPERATURE))?;
    let water =
        SaturatedWaterCavitation::try_new(pair, pair.surface_tension()?.surface_tension().value())?;
    println!(
        "case\ttime_s\tradius_m\tinterface_temperature_k\tinterface_velocity_m_s\tliquid_velocity_m_s\tvapor_mass_kg\tvapor_mass_rate_kg_s\tbubble_volume_rate_m3_s\tliquid_volume_flow_m3_s\tarrival_time_s\treceiver_pressure_pa\tcumulative_heat_j\tenergy_residual_j\tmass_residual_kg"
    );
    run("zero_heat", &water, 0.0)?;
    run("heat_pulse", &water, 0.001)?;
    Ok(())
}

fn run(
    name: &str,
    water: &SaturatedWaterCavitation,
    peak_heat_w: f64,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut bubble = MaterialVaporBubble::new(
        water,
        VaporBubbleState::try_new(0.0, RADIUS, 0.0, TEMPERATURE)?,
        VaporBubbleLimits {
            bubble: BubbleLimits {
                minimum_radius_m: 0.6 * RADIUS,
                maximum_radius_m: 1.5 * RADIUS,
                maximum_wall_mach_number: 0.01,
                maximum_vapor_density_ratio: 0.01,
            },
            minimum_temperature_k: TEMPERATURE - 2.0,
            maximum_temperature_k: TEMPERATURE + 2.0,
            maximum_pressure_ratio: 0.001,
            maximum_vapor_mach_number: 0.01,
            maximum_frozen_liquid_relative_change: 0.02,
            maximum_neglected_enthalpy_fraction: 0.001,
        },
    )?;
    let observer = bubble.monopole_observer(
        0.1,
        100_000.0,
        MonopoleLimits {
            maximum_ka: 0.1,
            minimum_distance_radii: 10.0,
            maximum_pressure_ratio: 0.01,
        },
    )?;
    let numerics = VaporBubbleNumerics {
        radial: BubbleNumerics {
            relative_tolerance: 1e-9,
            radius_absolute_tolerance_m: 1e-14,
            velocity_absolute_tolerance_m_s: 1e-9,
            initial_step_s: 2e-8,
            minimum_step_s: 1e-13,
            maximum_step_s: 2e-7,
            maximum_attempts: 20_000,
        },
        temperature_absolute_tolerance_k: 1e-9,
        energy_absolute_tolerance_j: 1e-16,
        mass_absolute_tolerance_kg: 1e-22,
        energy_balance_tolerance_j: 1e-12,
        mass_balance_tolerance_kg: 1e-19,
    };
    let equilibrium =
        water.state().saturation_pressure_pa() - 2.0 * water.state().surface_tension_n_m() / RADIUS;
    let pressure =
        |time: f64| equilibrium - 1000.0 * (core::f64::consts::TAU * time / DURATION).sin();
    let heat = |time: f64, _: &VaporBubbleThermalState| -> Result<f64, VaporBubbleError> {
        Ok(peak_heat_w * (core::f64::consts::PI * time / DURATION).sin().powi(2))
    };
    let initial = bubble.observe(equilibrium, 0.0)?;
    let initial_energy = initial.vapor_mass_kg * initial.vapor.vapor_internal_energy_j_kg();
    let mut total_heat = 0.0;
    let mut total_enthalpy = 0.0;
    let mut total_work = 0.0;
    let mut total_mass = 0.0;
    let mut maximum_energy_residual: f64 = 0.0;
    let mut maximum_mass_residual: f64 = 0.0;
    let mut maximum_pressure: f64 = 0.0;
    for i in 1..=SAMPLES {
        let step = bubble.advance_observed(
            i as f64 * DURATION / SAMPLES as f64,
            &numerics,
            &pressure,
            &heat,
            &observer,
            &mut || false,
        )?;
        let state = step.bubble.observation;
        total_heat += step.bubble.supplied_heat_j;
        total_enthalpy += step.bubble.transferred_enthalpy_j;
        total_work += step.bubble.pressure_work_j;
        total_mass += step.bubble.transferred_vapor_mass_kg;
        // Independent cumulative original-source balances across all windows.
        let energy_residual = state.vapor_mass_kg * state.vapor.vapor_internal_energy_j_kg()
            - initial_energy
            - total_heat
            - total_enthalpy
            + total_work;
        let mass_residual = state.vapor_mass_kg - initial.vapor_mass_kg - total_mass;
        maximum_energy_residual = maximum_energy_residual.max(energy_residual.abs());
        maximum_mass_residual = maximum_mass_residual.max(mass_residual.abs());
        maximum_pressure = maximum_pressure.max(step.acoustic.pressure_pa.abs());
        println!(
            "{name}\t{:.12e}\t{:.12e}\t{:.12e}\t{:.12e}\t{:.12e}\t{:.12e}\t{:.12e}\t{:.12e}\t{:.12e}\t{:.12e}\t{:.12e}\t{:.12e}\t{:.12e}\t{:.12e}",
            state.state.time_s(),
            state.state.radius_m(),
            state.state.temperature_k(),
            state.wall_velocity_m_s,
            state.state.liquid_velocity_m_s(),
            state.vapor_mass_kg,
            state.vapor_mass_rate_kg_s,
            state.volume_rate_m3_s,
            state.liquid_volume_flow_m3_s,
            step.acoustic.arrival_time_s,
            step.acoustic.pressure_pa,
            total_heat,
            energy_residual,
            mass_residual
        );
    }
    eprintln!(
        "{name}: peak_heat_w={peak_heat_w:e}; total_heat_j={total_heat:e}; peak_receiver_pressure_pa={maximum_pressure:e}; maximum_cumulative_energy_residual_j={maximum_energy_residual:e}; maximum_cumulative_mass_residual_kg={maximum_mass_residual:e}"
    );
    Ok(())
}
