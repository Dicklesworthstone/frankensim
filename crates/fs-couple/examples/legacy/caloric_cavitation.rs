//! Nonlinear gas heating/cooling in a cavitating bubble and its emitted sound.
//!
//! Run from the repository root:
//! cargo run -p fs-couple --features caloric-cavitation --example caloric_cavitation --
//! crates/fs-material/examples/material_gas_spring/mixture-a.tsv 2500 100
//!
//! GAS.tsv is the existing seven-column SYNTHETIC mixture example format, not
//! measured gas data. The same finite gas inventory follows two trajectories:
//! adiabatic and a prescribed heat-transfer coefficient h [W/(m² K)] with
//! Q = 4 pi R² h (T_liquid - T_gas). The heat law is an authored lumped closure,
//! not a resolved gas/liquid thermal boundary layer or a calibrated prediction.
//! Temperature-dependent NASA-9 calorics and transport remain source-owned.
//! IF97 supplies liquid water at 293.15 K; surface tension is authored 0.072 N/m.
//! The liquid is an isothermal reservoir, and vapor mass/latent heat are omitted
//! under explicit dilute-vapor admission. The bubble equation already includes
//! first-order acoustic back-reaction; the receiver removes no additional power.

#[path = "../../../fs-material/examples/material_gas_spring/source.rs"]
mod source;

use fs_aeroac::monopole::MonopoleLimits;
use fs_couple::cavitation::bubble::caloric::MaterialCaloricBubble;
use fs_flux::cavitation::caloric::{
    CaloricBubbleEquation, CaloricBubbleError, CaloricBubbleLimits, CaloricBubbleNumerics,
    CaloricGasObservation,
};
use fs_flux::cavitation::{
    BubbleLimits, BubbleNumerics, BubblePressureSample, BubbleState, KellerMiksisLimits,
};
use fs_material::cavitation::SaturatedWaterCavitation;
use fs_material::thermochemical_gas::energy::{GasEnergySolve, parcel::ClosedGasParcel};
use fs_qty::Temperature;
use fs_thermochem::water::if97::If97SaturationPair;
use std::io::Read;

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
const USAGE: &str =
    "usage: caloric_cavitation GAS.tsv [heat_transfer_W_m2_K=2500] [output_steps=100]";
const INITIAL_T: f64 = 293.15;
const RADIUS: f64 = 100e-6;
const AMBIENT: f64 = 101_325.0;
const DURATION: f64 = 500e-6;

fn main() -> Result<()> {
    let arguments: Vec<_> = std::env::args().skip(1).take(4).collect();
    if arguments == ["--help"] {
        println!("{USAGE}");
        return Ok(());
    }
    if arguments.is_empty() || arguments.len() > 3 {
        return Err(USAGE.into());
    }
    let transfer: f64 = arguments
        .get(1)
        .map(|s| s.parse())
        .transpose()?
        .unwrap_or(2500.0);
    let steps: u32 = arguments
        .get(2)
        .map(|s| s.parse())
        .transpose()?
        .unwrap_or(100);
    if !transfer.is_finite() || transfer <= 0.0 || !(1..=10_000).contains(&steps) {
        return Err(
            "positive finite heat-transfer coefficient and 1..=10000 output steps required".into(),
        );
    }
    let mut text = String::new();
    std::fs::File::open(&arguments[0])?
        .take(65_537)
        .read_to_string(&mut text)?;
    let gas_model = source::load(&text, &arguments[0])?;
    let water = SaturatedWaterCavitation::try_new(
        If97SaturationPair::at_temperature(Temperature::new(INITIAL_T))?,
        0.072,
    )?;
    let volume = (4.0 / 3.0) * std::f64::consts::PI * RADIUS.powi(3);
    let gas_pressure = AMBIENT - water.state().saturation_pressure_pa()
        + 2.0 * water.state().surface_tension_n_m() / RADIUS;
    let gas = ClosedGasParcel::new(
        gas_model,
        source::point(INITIAL_T, gas_pressure)?,
        volume,
    )?;
    let caloric = GasEnergySolve {
        minimum_temperature_k: 200.0,
        maximum_temperature_k: 1800.0,
        energy_tolerance_j_kg: 1e-6,
        maximum_iterations: 100,
    };
    let equation = CaloricBubbleEquation::KellerMiksis(KellerMiksisLimits {
        bubble: BubbleLimits {
            minimum_radius_m: 0.6 * RADIUS,
            maximum_radius_m: 1.4 * RADIUS,
            maximum_wall_mach_number: 0.01,
            maximum_vapor_density_ratio: 0.01,
        },
        maximum_pressure_ratio: 0.01,
    });
    let limits = CaloricBubbleLimits {
        maximum_gas_liquid_density_ratio: 0.01,
        maximum_vapor_pressure_fraction: 0.1,
        maximum_gas_wall_mach_number: 0.01,
    };
    let omega = std::f64::consts::TAU / DURATION;
    // A smooth compression/release cycle. Its exact derivative is part of KM.
    let pressure = |time_s: f64| BubblePressureSample {
        absolute_pressure_pa: AMBIENT + 10_000.0 * (1.0 - (omega * time_s).cos()),
        time_derivative_pa_s: 10_000.0 * omega * (omega * time_s).sin(),
    };
    let initial_state = BubbleState::try_new(0.0, RADIUS, 0.0)?;
    let make_bubble = || {
        MaterialCaloricBubble::new(
            &water,
            gas.clone(),
            equation,
            initial_state,
            pressure(0.0),
            limits,
            caloric,
        )
    };
    let mut adiabatic = make_bubble()?;
    let mut exchanging = make_bubble()?;
    let observer = adiabatic.monopole_observer(
        0.1,
        100_000.0,
        MonopoleLimits {
            maximum_ka: 0.1,
            minimum_distance_radii: 10.0,
            maximum_pressure_ratio: 0.01,
        },
    )?;
    let numerics = CaloricBubbleNumerics {
        radial: BubbleNumerics {
            relative_tolerance: 1e-8,
            radius_absolute_tolerance_m: 1e-13,
            velocity_absolute_tolerance_m_s: 1e-8,
            initial_step_s: 1e-8,
            minimum_step_s: 1e-12,
            maximum_step_s: 2.5e-7,
            maximum_attempts: 50_000,
        },
        energy_absolute_tolerance_j: 1e-17,
        energy_balance_tolerance_j: 1e-15,
    };
    let no_heat = |_: f64, _: &CaloricGasObservation| Ok::<_, CaloricBubbleError>(0.0);
    let heat = |_: f64, state: &CaloricGasObservation| {
        Ok::<_, CaloricBubbleError>(
            4.0 * std::f64::consts::PI * state.radius_m.powi(2)
                * transfer
                * (INITIAL_T - state.temperature_k),
        )
    };
    let mut heat_totals = [0.0; 2];
    let mut work_totals = [0.0; 2];
    println!("# IF97 water, synthetic caller-supplied gas; isothermal liquid reservoir");
    println!(
        "# radius=100 um, receiver=0.1 m, compression pulse=20 kPa, h={transfer} W/(m2 K)"
    );
    println!(
        "model,time_us,arrival_us,radius_um,gas_temperature_K,gas_pressure_Pa,receiver_pressure_Pa,heat_J,gas_work_J,gas_energy_residual_J"
    );
    for index in 1..=steps {
        let end_time = DURATION * f64::from(index) / f64::from(steps);
        let outputs = [
            adiabatic.advance_observed(
                end_time, &numerics, &pressure, &no_heat, &observer, &mut || false,
            )?,
            exchanging.advance_observed(
                end_time, &numerics, &pressure, &heat, &observer, &mut || false,
            )?,
        ];
        for (model_index, (label, output)) in
            ["adiabatic", "heat_exchange"].into_iter().zip(outputs).enumerate()
        {
            let state = output.bubble.observation.gas;
            heat_totals[model_index] += output.bubble.supplied_heat_j;
            work_totals[model_index] += output.bubble.work_by_gas_j;
            println!(
                "{label},{:.6},{:.6},{:.9},{:.9},{:.9},{:.9e},{:.9e},{:.9e},{:.3e}",
                state.time_s * 1e6,
                output.acoustic.arrival_time_s * 1e6,
                state.radius_m * 1e6,
                state.temperature_k,
                state.gas_pressure_pa,
                output.acoustic.pressure_pa,
                heat_totals[model_index],
                work_totals[model_index],
                output.bubble.gas_energy_residual_j,
            );
        }
    }
    for (index, (label, owner)) in ["adiabatic", "heat_exchange"]
        .into_iter()
        .zip([adiabatic, exchanging])
        .enumerate()
    {
        let change = gas.stored_energy_change_to(owner.gas(), || false)?;
        println!(
            "# {label}: original-law gas energy change={change:.9e} J; accumulated balance residual={:.3e} J",
            change - heat_totals[index] + work_totals[index],
        );
    }
    Ok(())
}
