//! Resolved liquid diffusion, evaporation/condensation and water-borne sound.
//!
//! Run: cargo run -p fs-couple --features cavitation --example conducting_vapor_cavitation
//!
//! A 100 um vapor bubble starts in 373.15 K water. A smooth 40 us pressure
//! cycle drives growth and collapse. Heat is computed from the retained liquid
//! field with original IF97 calorics/transport; no interface heat coefficient
//! or prescribed bubble heat pulse is supplied. The finite grid exchanges heat
//! and advected enthalpy with an explicit outer bath. Three graded grids and a
//! halved maximum step expose spatial and temporal sensitivity of this case.
//! This is a bounded spherical local-equilibrium model, not shock collapse or
//! experimental validation. The outer-domain and material admission limits
//! remain explicit and must be reconsidered for a different physical case.

use fs_aeroac::monopole::MonopoleLimits;
use fs_couple::cavitation::bubble::vapor::conducting::MaterialConductingVaporBubble;
use fs_flux::cavitation::vapor::conducting::ConductingVaporBubbleNumerics;
use fs_flux::cavitation::vapor::liquid::{
    LiquidThermalLimits, SphericalLiquidField, SphericalLiquidGrid,
};
use fs_flux::cavitation::vapor::{VaporBubbleLimits, VaporBubbleNumerics, VaporBubbleState};
use fs_flux::cavitation::{BubbleLimits, BubbleNumerics};
use fs_material::cavitation::SaturatedWaterCavitation;
use fs_qty::Temperature;
use fs_thermochem::water::if97::If97SaturationPair;

const RADIUS: f64 = 100e-6;
const TEMPERATURE: f64 = 373.15;
const DURATION: f64 = 40e-6;
const SAMPLES: usize = 160;
const LIQUID_EXTENT: f64 = 40e-6;

#[derive(Clone, Copy)]
struct Summary {
    peak_pressure_pa: f64,
    peak_heat_w: f64,
    final_radius_m: f64,
    final_vapor_mass_kg: f64,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let pair = If97SaturationPair::at_temperature(Temperature::new(TEMPERATURE))?;
    let water =
        SaturatedWaterCavitation::try_new(pair, pair.surface_tension()?.surface_tension().value())?;
    println!(
        "case\ttime_s\tradius_m\tinterface_temperature_k\tfirst_liquid_temperature_k\tinterface_velocity_m_s\tliquid_velocity_m_s\tvapor_mass_kg\tvapor_mass_rate_kg_s\tinterface_heat_w\tarrival_time_s\treceiver_pressure_pa\tcumulative_interface_heat_j\tcumulative_viscous_heat_j\tcumulative_energy_residual_j\tcumulative_mass_residual_kg"
    );
    let coarse = run("liquid_16", &water, 16, 2e-7)?;
    let medium = run("liquid_32", &water, 32, 2e-7)?;
    let fine = run("liquid_64", &water, 64, 2e-7)?;
    let time_fine = run("liquid_64_half_step", &water, 64, 1e-7)?;
    for (name, sample) in [
        ("16_vs_64", coarse),
        ("32_vs_64", medium),
        ("half_step_vs_64", time_fine),
    ] {
        eprintln!(
            "{name}: peak_pressure_relative_change={:e}; peak_heat_relative_change={:e}; final_radius_relative_change={:e}; final_vapor_mass_relative_change={:e}",
            relative(sample.peak_pressure_pa, fine.peak_pressure_pa),
            relative(sample.peak_heat_w, fine.peak_heat_w),
            relative(sample.final_radius_m, fine.final_radius_m),
            relative(sample.final_vapor_mass_kg, fine.final_vapor_mass_kg),
        );
    }
    Ok(())
}

fn relative(value: f64, reference: f64) -> f64 {
    (value - reference).abs() / reference.abs().max(1e-30)
}

fn run(
    name: &str,
    water: &SaturatedWaterCavitation,
    cells: usize,
    maximum_step_s: f64,
) -> Result<Summary, Box<dyn std::error::Error>> {
    let mut faces = [0.0; 65];
    let grading = 3.0_f64;
    for (i, face) in faces.iter_mut().enumerate().take(cells + 1) {
        *face = RADIUS
            + LIQUID_EXTENT * (grading * i as f64 / cells as f64).exp_m1() / grading.exp_m1();
    }
    let liquid = SphericalLiquidField::new(
        *water,
        SphericalLiquidGrid::from_face_radii(&faces[..=cells])?,
        LiquidThermalLimits {
            minimum_temperature_k: TEMPERATURE - 2.0,
            maximum_temperature_k: TEMPERATURE + 2.0,
            maximum_source_relative_change: 0.02,
            maximum_pressure_path_capacity_fraction: 0.01,
            maximum_pressure_heating_temperature_k: 0.01,
            maximum_outflow_peclet: 0.5,
        },
    )?;
    let anchor = liquid.reference_enthalpy_j_kg();
    let mut bubble = MaterialConductingVaporBubble::new(
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
        liquid,
        &[TEMPERATURE; 64][..cells],
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
    let numerics = ConductingVaporBubbleNumerics {
        bubble: VaporBubbleNumerics {
            radial: BubbleNumerics {
                relative_tolerance: 1e-9,
                radius_absolute_tolerance_m: 1e-14,
                velocity_absolute_tolerance_m_s: 1e-9,
                initial_step_s: 2e-8,
                minimum_step_s: 1e-13,
                maximum_step_s,
                maximum_attempts: 20_000,
            },
            temperature_absolute_tolerance_k: 1e-9,
            energy_absolute_tolerance_j: 1e-16,
            mass_absolute_tolerance_kg: 1e-22,
            energy_balance_tolerance_j: 1e-13,
            mass_balance_tolerance_kg: 1e-19,
        },
        liquid_temperature_absolute_tolerance_k: 1e-9,
        liquid_energy_balance_tolerance_j: 1e-13,
        total_energy_balance_tolerance_j: 1e-13,
    };
    let equilibrium =
        water.state().saturation_pressure_pa() - 2.0 * water.state().surface_tension_n_m() / RADIUS;
    let pressure =
        |time: f64| equilibrium - 1000.0 * (core::f64::consts::TAU * time / DURATION).sin();
    let initial = bubble.observe(equilibrium, TEMPERATURE)?;
    let mut total_heat = 0.0;
    let mut total_work = 0.0;
    let mut total_mass = 0.0;
    let mut total_outer_heat = 0.0;
    let mut total_outer_enthalpy = 0.0;
    let mut total_viscous = 0.0;
    let mut maximum_energy_residual: f64 = 0.0;
    let mut maximum_mass_residual: f64 = 0.0;
    let mut maximum_pressure: f64 = 0.0;
    let mut maximum_heat: f64 = 0.0;
    let mut accepted_steps = 0;
    let mut rhs_evaluations = 0;
    let mut last_mass = initial.bubble.vapor_mass_kg;
    for i in 1..=SAMPLES {
        let step = bubble.advance_observed(
            i as f64 * DURATION / SAMPLES as f64,
            &numerics,
            &pressure,
            &|_| TEMPERATURE,
            &observer,
            &mut || false,
        )?;
        let field = step.bubble;
        let phase = field.bubble;
        let state = phase.observation;
        total_heat += phase.supplied_heat_j;
        total_work += phase.pressure_work_j;
        total_mass += phase.transferred_vapor_mass_kg;
        total_outer_heat += field.outer_conductive_heat_j;
        total_outer_enthalpy += field.outer_advected_enthalpy_j;
        total_viscous += field.viscous_heat_j;
        accepted_steps += phase.accepted_steps;
        rhs_evaluations += phase.rhs_evaluations;
        // Independent cumulative balance, evaluated from original endpoint
        // calorics and every outer port, not the sum of per-window residuals.
        let phase_change = (state.vapor_mass_kg - initial.bubble.vapor_mass_kg)
            * (state.vapor.vapor_internal_energy_j_kg() - anchor)
            + initial.bubble.vapor_mass_kg
                * (state.vapor.vapor_internal_energy_j_kg()
                    - initial.bubble.vapor.vapor_internal_energy_j_kg());
        let liquid_change = field.liquid.energy_change_from(&initial.liquid)?;
        let energy_residual =
            phase_change + liquid_change + total_work + total_outer_heat + total_outer_enthalpy
                - total_viscous;
        let mass_residual = state.vapor_mass_kg - initial.bubble.vapor_mass_kg - total_mass;
        maximum_energy_residual = maximum_energy_residual.max(energy_residual.abs());
        maximum_mass_residual = maximum_mass_residual.max(mass_residual.abs());
        maximum_pressure = maximum_pressure.max(step.acoustic.pressure_pa.abs());
        maximum_heat = maximum_heat.max(field.liquid.interface_heat_power_w.abs());
        last_mass = state.vapor_mass_kg;
        println!(
            "{name}\t{:.12e}\t{:.12e}\t{:.12e}\t{:.12e}\t{:.12e}\t{:.12e}\t{:.12e}\t{:.12e}\t{:.12e}\t{:.12e}\t{:.12e}\t{:.12e}\t{:.12e}\t{:.12e}\t{:.12e}",
            state.state.time_s(),
            state.state.radius_m(),
            state.state.temperature_k(),
            bubble.liquid_temperatures_k()[0],
            state.wall_velocity_m_s,
            state.state.liquid_velocity_m_s(),
            state.vapor_mass_kg,
            state.vapor_mass_rate_kg_s,
            field.liquid.interface_heat_power_w,
            step.acoustic.arrival_time_s,
            step.acoustic.pressure_pa,
            total_heat,
            total_viscous,
            energy_residual,
            mass_residual,
        );
    }
    eprintln!(
        "{name}: cells={cells}; maximum_step_s={maximum_step_s:e}; peak_receiver_pressure_pa={maximum_pressure:e}; peak_interface_heat_w={maximum_heat:e}; total_interface_heat_j={total_heat:e}; maximum_cumulative_energy_residual_j={maximum_energy_residual:e}; maximum_cumulative_mass_residual_kg={maximum_mass_residual:e}; accepted_steps={accepted_steps}; rhs_evaluations={rhs_evaluations}"
    );
    Ok(Summary {
        peak_pressure_pa: maximum_pressure,
        peak_heat_w: maximum_heat,
        final_radius_m: bubble.state().radius_m(),
        final_vapor_mass_kg: last_mass,
    })
}
