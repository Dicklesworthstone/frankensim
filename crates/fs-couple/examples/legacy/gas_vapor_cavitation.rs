//! A conserved gas nucleus with resolved liquid heat, water phase exchange and sound.
//!
//! Run:
//! cargo run -p fs-couple --features caloric-cavitation --example gas_vapor_cavitation -- \
//! crates/fs-material/examples/material_gas_spring/mixture-a.tsv 0.02 2e-5
//!
//! GAS.tsv uses the existing SYNTHETIC source-bearing example format. The two
//! required numeric arguments explicitly author lower bounds on WET-mixture
//! conductivity [W/(m K)] and water/gas binary diffusivity [m2/s] over the
//! printed temperature, pressure and composition domain. They are not inferred
//! from dry-gas transport, fitted to sound, or presented as measured data.
//!
//! A 10 um nucleus in 353.15 K water starts at mechanical equilibrium. A
//! smooth rarefaction/compression cycle changes water vapor mass while every
//! noncondensable species retains its original inventory. One retained liquid
//! field supplies heat; the one native trajectory supplies the acoustic source.
//! This example exercises a bounded common-temperature, rapid-mixing model,
//! not spatial gas diffusion, finite-rate interface kinetics or shock collapse.

#[path = "../../../fs-material/examples/material_gas_spring/source.rs"]
mod source;

use fs_aeroac::monopole::MonopoleLimits;
use fs_couple::cavitation::bubble::vapor::conducting::MaterialConductingVaporBubble;
use fs_flux::cavitation::vapor::conducting::ConductingVaporBubbleNumerics;
use fs_flux::cavitation::vapor::conducting::gas::{GasVaporLimits, GasVaporTransportBounds};
use fs_flux::cavitation::vapor::liquid::{
    LiquidThermalLimits, SphericalLiquidField, SphericalLiquidGrid,
};
use fs_flux::cavitation::vapor::{VaporBubbleLimits, VaporBubbleNumerics, VaporBubbleState};
use fs_flux::cavitation::{BubbleLimits, BubbleNumerics};
use fs_material::cavitation::SaturatedWaterCavitation;
use fs_material::cavitation::vapor::gas::ConservedBubbleGas;
use fs_material::thermochemical_gas::energy::{GasEnergySolve, parcel::ClosedGasParcel};
use fs_qty::Temperature;
use fs_thermochem::water::if97::If97SaturationPair;
use std::io::Read;

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
const USAGE: &str =
    "usage: gas_vapor_cavitation GAS.tsv MIXTURE_K_MIN_W_m_K MIXTURE_D_MIN_m2_s [output_steps=200]";
const RADIUS: f64 = 10e-6;
const TEMPERATURE: f64 = 353.15;
const AMBIENT: f64 = 101_325.0;
const AMPLITUDE: f64 = 40_000.0;
const DURATION: f64 = 4e-3;
const CELLS: usize = 16;

fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).take(5).collect();
    if args == ["--help"] {
        println!("{USAGE}");
        return Ok(());
    }
    if args.len() < 3 || args.len() > 4 {
        return Err(USAGE.into());
    }
    let conductivity: f64 = args[1].parse()?;
    let diffusivity: f64 = args[2].parse()?;
    let samples: u32 = args.get(3).map(|s| s.parse()).transpose()?.unwrap_or(200);
    if !conductivity.is_finite()
        || conductivity <= 0.0
        || !diffusivity.is_finite()
        || diffusivity <= 0.0
        || !(1..=10_000).contains(&samples)
    {
        return Err(
            "positive finite wet-mixture bounds and 1..=10000 output steps required".into(),
        );
    }
    let mut text = String::new();
    std::fs::File::open(&args[0])?
        .take(65_537)
        .read_to_string(&mut text)?;
    let gas_model = source::load(&text, &args[0])?;
    let pair = If97SaturationPair::at_temperature(Temperature::new(TEMPERATURE))?;
    let water =
        SaturatedWaterCavitation::try_new(pair, pair.surface_tension()?.surface_tension().value())?;
    let initial_volume = (4.0 / 3.0) * core::f64::consts::PI * RADIUS.powi(3);
    let initial_gas_pressure = AMBIENT - water.state().saturation_pressure_pa()
        + 2.0 * water.state().surface_tension_n_m() / RADIUS;
    let gas = ConservedBubbleGas::new(
        ClosedGasParcel::new(
            gas_model,
            source::point(TEMPERATURE, initial_gas_pressure)?,
            initial_volume,
        )?,
        GasEnergySolve {
            minimum_temperature_k: TEMPERATURE - 2.0,
            maximum_temperature_k: TEMPERATURE + 2.0,
            energy_tolerance_j_kg: 1e-7,
            maximum_iterations: 100,
        },
    )?;
    let mut faces = [0.0; CELLS + 1];
    for (i, face) in faces.iter_mut().enumerate() {
        *face = RADIUS * (1.0 + 2.0 * i as f64 / CELLS as f64);
    }
    let field = SphericalLiquidField::new(
        water,
        SphericalLiquidGrid::from_face_radii(&faces)?,
        LiquidThermalLimits {
            minimum_temperature_k: TEMPERATURE - 2.0,
            maximum_temperature_k: TEMPERATURE + 2.0,
            maximum_source_relative_change: 0.05,
            maximum_pressure_path_capacity_fraction: 0.01,
            maximum_pressure_heating_temperature_k: 0.01,
            maximum_outflow_peclet: 0.5,
        },
    )?;
    let anchor = field.reference_enthalpy_j_kg();
    let omega = core::f64::consts::TAU / DURATION;
    let gas_limits = GasVaporLimits {
        transport: GasVaporTransportBounds {
            minimum_temperature_k: TEMPERATURE - 2.0,
            maximum_temperature_k: TEMPERATURE + 2.0,
            minimum_total_pressure_pa: 50_000.0,
            maximum_total_pressure_pa: 200_000.0,
            minimum_vapor_mass_fraction: 0.01,
            maximum_vapor_mass_fraction: 0.95,
            minimum_thermal_conductivity_w_m_k: conductivity,
            minimum_binary_diffusivity_m2_s: diffusivity,
        },
        maximum_total_density_ratio: 0.01,
        maximum_mixture_mach_number: 0.01,
        maximum_pressure_relaxation_ratio: 0.01,
        maximum_internal_temperature_difference_k: 1.0,
        maximum_relative_composition_gradient: 0.05,
        maximum_thermal_response_ratio: 0.1,
        maximum_diffusion_response_ratio: 0.1,
        // sin^3 has a third harmonic. This bound covers its imposed bandwidth
        // as well as the maximum relative pressure forcing rate.
        maximum_imposed_rate_per_s: 3.0 * omega,
    };
    let mut bubble = MaterialConductingVaporBubble::new(
        &water,
        VaporBubbleState::try_new(0.0, RADIUS, 0.0, TEMPERATURE)?,
        VaporBubbleLimits {
            bubble: BubbleLimits {
                minimum_radius_m: 0.6 * RADIUS,
                maximum_radius_m: 1.8 * RADIUS,
                maximum_wall_mach_number: 0.01,
                maximum_vapor_density_ratio: 0.01,
            },
            minimum_temperature_k: TEMPERATURE - 2.0,
            maximum_temperature_k: TEMPERATURE + 2.0,
            maximum_pressure_ratio: 0.001,
            maximum_vapor_mach_number: 0.01,
            maximum_frozen_liquid_relative_change: 0.05,
            maximum_neglected_enthalpy_fraction: 0.005,
        },
        field,
        &[TEMPERATURE; CELLS],
    )?
    .with_noncondensable_gas(gas.clone(), gas_limits)?;
    let observer = bubble.monopole_observer(
        0.1,
        100_000.0,
        MonopoleLimits {
            maximum_ka: 0.1,
            minimum_distance_radii: 10.0,
            maximum_pressure_ratio: 0.01,
        },
    )?;
    let controls = ConductingVaporBubbleNumerics {
        bubble: VaporBubbleNumerics {
            radial: BubbleNumerics {
                relative_tolerance: 1e-8,
                radius_absolute_tolerance_m: 1e-15,
                velocity_absolute_tolerance_m_s: 1e-8,
                initial_step_s: 1e-8,
                minimum_step_s: 1e-14,
                maximum_step_s: 1e-6,
                maximum_attempts: 50_000,
            },
            temperature_absolute_tolerance_k: 1e-8,
            energy_absolute_tolerance_j: 1e-19,
            mass_absolute_tolerance_kg: 1e-24,
            energy_balance_tolerance_j: 1e-15,
            mass_balance_tolerance_kg: 1e-21,
        },
        liquid_temperature_absolute_tolerance_k: 1e-8,
        liquid_energy_balance_tolerance_j: 1e-15,
        total_energy_balance_tolerance_j: 1e-15,
    };
    let pressure = |time: f64| AMBIENT - AMPLITUDE * (omega * time).sin().powi(3);
    let initial = bubble.observe(AMBIENT, TEMPERATURE)?;
    let mut pressure_work = 0.0;
    let mut outer_heat = 0.0;
    let mut outer_enthalpy = 0.0;
    let mut viscous_heat = 0.0;
    let mut vapor_transfer = 0.0;
    let mut peak_pressure: f64 = 0.0;
    let mut maximum_energy_residual: f64 = 0.0;
    let mut maximum_mass_residual: f64 = 0.0;
    let mut maximum_gas_energy_residual: f64 = 0.0;
    let mut minimum_radius = RADIUS;
    let mut maximum_radius = RADIUS;
    let mut minimum_vapor_mass = initial.bubble.vapor_mass_kg;
    let mut maximum_vapor_mass = initial.bubble.vapor_mass_kg;
    let mut gas_energy_integral = 0.0;
    let mut evaporation_samples = 0;
    let mut condensation_samples = 0;
    eprintln!(
        "SYNTHETIC gas; authored wet-mixture lower bounds k={conductivity:e} W/(m K), D={diffusivity:e} m2/s; T={}..{} K; total p=50000..200000 Pa; vapor mass fraction=0.01..0.95; fixed gas mass={:e} kg",
        TEMPERATURE - 2.0,
        TEMPERATURE + 2.0,
        gas.mass_kg(),
    );
    println!(
        "time_s,radius_m,temperature_k,vapor_pressure_pa,gas_partial_pressure_pa,total_bubble_pressure_pa,gas_mass_kg,vapor_mass_kg,vapor_mass_rate_kg_s,interface_heat_w,arrival_time_s,receiver_pressure_pa,gas_energy_change_j,gas_work_j,thermal_energy_residual_j,vapor_mass_residual_kg,gas_energy_residual_j,internal_temperature_scale_k,species_gradient_scale"
    );
    let mut gas_work = 0.0;
    for i in 1..=samples {
        let end = DURATION * f64::from(i) / f64::from(samples);
        let out = bubble.advance_observed(
            end,
            &controls,
            &pressure,
            &|_| TEMPERATURE,
            &observer,
            &mut || false,
        )?;
        let step = out.bubble;
        let state = step.bubble.observation;
        let gas_state = step
            .gas
            .ok_or("mixed endpoint lost the original gas inventory")?;
        if gas_state.mass_kg.to_bits() != gas.mass_kg().to_bits() {
            return Err("noncondensable inventory changed".into());
        }
        pressure_work += step.bubble.pressure_work_j;
        outer_heat += step.outer_conductive_heat_j;
        outer_enthalpy += step.outer_advected_enthalpy_j;
        viscous_heat += step.viscous_heat_j;
        vapor_transfer += step.bubble.transferred_vapor_mass_kg;
        gas_energy_integral += step.integrated_noncondensable_energy_j;
        gas_work += step.noncondensable_pressure_work_j;
        let phase_energy_change = (state.vapor_mass_kg - initial.bubble.vapor_mass_kg)
            * (state.vapor.vapor_internal_energy_j_kg() - anchor)
            + initial.bubble.vapor_mass_kg
                * (state.vapor.vapor_internal_energy_j_kg()
                    - initial.bubble.vapor.vapor_internal_energy_j_kg());
        let gas_energy_change = gas.energy_change_j(TEMPERATURE, state.state.temperature_k())?;
        let liquid_change = step.liquid.energy_change_from(&initial.liquid)?;
        let energy_residual = phase_energy_change
            + gas_energy_change
            + liquid_change
            + pressure_work
            + outer_heat
            + outer_enthalpy
            - viscous_heat;
        let mass_residual = state.vapor_mass_kg - initial.bubble.vapor_mass_kg - vapor_transfer;
        let gas_energy_residual = gas_energy_change - gas_energy_integral;
        maximum_energy_residual = maximum_energy_residual.max(energy_residual.abs());
        maximum_mass_residual = maximum_mass_residual.max(mass_residual.abs());
        maximum_gas_energy_residual = maximum_gas_energy_residual.max(gas_energy_residual.abs());
        peak_pressure = peak_pressure.max(out.acoustic.pressure_pa.abs());
        minimum_radius = minimum_radius.min(state.state.radius_m());
        maximum_radius = maximum_radius.max(state.state.radius_m());
        minimum_vapor_mass = minimum_vapor_mass.min(state.vapor_mass_kg);
        maximum_vapor_mass = maximum_vapor_mass.max(state.vapor_mass_kg);
        evaporation_samples += usize::from(state.vapor_mass_rate_kg_s > 0.0);
        condensation_samples += usize::from(state.vapor_mass_rate_kg_s < 0.0);
        println!(
            "{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e}",
            state.state.time_s(),
            state.state.radius_m(),
            state.state.temperature_k(),
            state.vapor.pressure_pa(),
            state.noncondensable_pressure_pa,
            state.bubble_pressure_pa,
            gas_state.mass_kg,
            state.vapor_mass_kg,
            state.vapor_mass_rate_kg_s,
            step.liquid.interface_heat_power_w,
            out.acoustic.arrival_time_s,
            out.acoustic.pressure_pa,
            gas_energy_change,
            gas_work,
            energy_residual,
            mass_residual,
            gas_energy_residual,
            gas_state.estimated_internal_temperature_difference_k,
            gas_state.relative_composition_gradient_scale,
        );
    }
    eprintln!(
        "completed_samples={samples}; evaporation_samples={evaporation_samples}; condensation_samples={condensation_samples}; radius_range_m={minimum_radius:e}..{maximum_radius:e}; vapor_mass_range_kg={minimum_vapor_mass:e}..{maximum_vapor_mass:e}; unchanged_gas_mass_kg={:e}; peak_sampled_receiver_pressure_pa={peak_pressure:e}; maximum_cumulative_thermal_energy_residual_j={maximum_energy_residual:e}; maximum_cumulative_vapor_mass_residual_kg={maximum_mass_residual:e}; maximum_cumulative_gas_energy_residual_j={maximum_gas_energy_residual:e}",
        gas.mass_kg(),
    );
    Ok(())
}
