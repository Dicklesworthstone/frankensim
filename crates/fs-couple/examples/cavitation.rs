//! Source-bound liquid-film rupture/refill and a resolved bubble sound source.
//!
//! Run with `cargo run -p fs-couple --features cavitation --example cavitation`.
//! Water saturation/transport comes from the original IF97 owner. The interface
//! tension is explicitly prescribed; temperature remains fixed throughout.
//! Film content is not vapor mass. Keller--Miksis supplies first-order liquid
//! compressibility; the compact observer adds no second radiation loss.

use fs_aeroac::monopole::MonopoleLimits;
use fs_couple::cavitation::bubble::MaterialBubble;
use fs_couple::cavitation::{
    CavitationFilmEnvironment, CavitationSolveConfig, MaterialCavitatingFilm,
};
use fs_flux::cavitation::{
    BubbleGas, BubbleLimits, BubbleNumerics, BubblePressureSample, BubbleState, KellerMiksisLimits,
};
use fs_material::cavitation::SaturatedWaterCavitation;
use fs_qty::Temperature;
use fs_thermochem::water::if97::If97SaturationPair;
use fs_tribo::resistive_film::{FilmCell, FilmChannel, FilmLimits, GapPort};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let water = SaturatedWaterCavitation::try_new(
        If97SaturationPair::at_temperature(Temperature::new(293.15))?,
        0.072,
    )?;
    let ambient_pa = 101_325.0;
    let gap = |reference_m, closure| GapPort {
        reference_m,
        closure: vec![closure],
    };
    let film = MaterialCavitatingFilm::new(
        &water,
        vec![FilmCell {
            area_m2: 1e-4,
            gap: gap(20e-6, 1.0),
        }],
        vec![FilmChannel {
            from: 0,
            to: None,
            width_m: 1e-3,
            length_m: 1e-2,
            gap: gap(10e-6, 0.0),
        }],
        1,
        FilmLimits {
            minimum_cell_gap_m: 1e-6,
            maximum_gap_m: 1e-3,
            maximum_pressure_pa: 2e5,
        },
        CavitationFilmEnvironment {
            ambient_pressure_pa: ambient_pa,
            ambient_liquid_fraction: 1.0,
            maximum_vapor_density_ratio: 0.01,
            solve: CavitationSolveConfig {
                maximum_iterations: 32,
                absolute_mass_tolerance_kg: 1e-17,
                relative_mass_tolerance: 1e-11,
            },
        },
    )?;
    let initial = film.initial_state(&[0.0])?;
    let opened = film.advance(&initial, &[-20e-6], &[0.0], 0.01, || false)?;
    let refilled = film.advance(&opened.state, &[0.0], &[0.0], 0.01, || false)?;
    println!("IF97 water at 293.15 K; prescribed surface tension 0.072 N/m");
    println!(
        "saturation pressure: {:.3} Pa",
        water.state().saturation_pressure_pa()
    );
    for (name, step) in [("opened", &opened), ("refilled", &refilled)] {
        println!(
            "film {name}: liquid content={:.6}, absolute pressure={:.3} Pa, mass residual={:.3e} kg",
            step.state.liquid_fraction()[0],
            step.state.pressure_pa()[0],
            step.report.mass_residual_kg,
        );
    }

    let radius_m = 100e-6;
    let gas = BubbleGas::equilibrium(water.state(), radius_m, ambient_pa, 1.4)?;
    let omega = 2.0 * std::f64::consts::PI * 2000.0;
    let pressure = |t: f64| BubblePressureSample {
        absolute_pressure_pa: ambient_pa + 200.0 * (omega * t).sin(),
        time_derivative_pa_s: 200.0 * omega * (omega * t).cos(),
    };
    let mut bubble = MaterialBubble::new_compressible(
        &water,
        gas,
        KellerMiksisLimits {
            bubble: BubbleLimits {
                minimum_radius_m: 0.5 * radius_m,
                maximum_radius_m: 2.0 * radius_m,
                maximum_wall_mach_number: 0.01,
                maximum_vapor_density_ratio: 0.01,
            },
            maximum_pressure_ratio: 0.01,
        },
        BubbleState::try_new(0.0, radius_m, 0.0)?,
        pressure(0.0),
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
    let numerics = BubbleNumerics {
        relative_tolerance: 1e-8,
        radius_absolute_tolerance_m: 1e-13,
        velocity_absolute_tolerance_m_s: 1e-8,
        initial_step_s: 1e-8,
        minimum_step_s: 1e-12,
        maximum_step_s: 1e-7,
        maximum_attempts: 10_000,
    };
    // A weak drive keeps this example in the declared spherical, weakly
    // compressible and linear-radiation domains. The analytic drive derivative
    // above is required by Keller--Miksis. Observation bandwidth includes the
    // natural small-amplitude response as well as this imposed frequency.
    println!("Keller--Miksis bubble with material-owned liquid compressibility");
    println!("source_us arrival_us radius_um acoustic_pressure_Pa (observer at 0.1 m)");
    for endpoint in [25e-6, 50e-6, 75e-6, 100e-6] {
        let step =
            bubble.advance_observed(endpoint, &numerics, &pressure, &observer, &mut || false)?;
        println!(
            "{:.3} {:.3} {:.6} {:.6e}",
            step.acoustic.source_time_s * 1e6,
            step.acoustic.arrival_time_s * 1e6,
            step.bubble.observation.state.radius_m() * 1e6,
            step.acoustic.pressure_pa,
        );
    }
    Ok(())
}
