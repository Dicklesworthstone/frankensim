//! Source-bound liquid-film rupture/refill and a resolved bubble sound source.
//!
//! Run with `cargo run -p fs-couple --features cavitation --example cavitation`.
//! Water saturation/transport comes from the original IF97 owner. The interface
//! tension is explicitly prescribed; temperature remains fixed throughout.
//! Film content is not vapor mass. Keller--Miksis supplies first-order liquid
//! compressibility; the compact observer adds no second radiation loss.
//! A final pressure sweep preserves nucleus gas inventories and number counts
//! while updating equilibrium radii and coherent acoustic propagation.
//! A gas/liquid diffusion spectrum then separates gas conduction, liquid
//! conduction and viscous absorption from scattering, using original IF97
//! water and dry-air transport at each nucleus's actual gas pressure.

use fs_aeroac::curle3d::moving::{MovingCompactLimits, MovingCompactObserver};
use fs_aeroac::monopole::MonopoleLimits;
use fs_couple::cavitation::bubble::MaterialBubble;
use fs_couple::cavitation::bubbly::thermal::{
    ConductingLiquidPopulationLimits, ThermalPolydisperseBubblyLiquid,
};
use fs_couple::cavitation::bubbly::{BubblyLiquidLimits, PolydisperseBubblyLiquid};
use fs_couple::cavitation::cluster::MaterialBubbleCluster;
use fs_couple::cavitation::{
    CavitationFilmEnvironment, CavitationSolveConfig, MaterialCavitatingFilm,
};
use fs_couple::pcm_wav::observation::moving::{
    CompactSourceHistory, MovingCompactRenderer, MovingPressureBudget, MovingPressureConfig,
};
use fs_flux::cavitation::inception::{
    EquilibriumBranch, EquilibriumControls, InceptionError, blake_threshold,
    equilibrium_at_pressure,
};
use fs_flux::cavitation::thermal::ThermalBubbleLimits;
use fs_flux::cavitation::{
    BubbleCluster, BubbleClusterLimits, BubbleClusterNumerics, BubbleGas, BubbleLimits,
    BubbleNumerics, BubblePopulation, BubblePopulationClass, BubblePressureSample, BubbleState,
    ClusterBubble, KellerMiksisBubble, KellerMiksisLimits, RayleighPlessetBubble,
};
use fs_material::cavitation::SaturatedWaterCavitation;
use fs_material::fluid_acoustics::FluidAcoustics;
use fs_material::gas::{GasSpec, GasState};
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
    cluster_noise(&water)?;
    nucleus_inception(&water)?;
    pressure_sweep(&water)
}

// Continue the SAME gas inventories and number densities through quasistatic
// operating pressures. Liquid coefficients remain frozen at the authored
// temperature; this is a held-volume distribution, not a closed-parcel process.
fn pressure_sweep(water: &SaturatedWaterCavitation) -> Result<(), Box<dyn std::error::Error>> {
    let liquid = water.state();
    let reference_pressure = 101_325.0;
    let mut classes = Vec::with_capacity(2);
    for (radius, fraction) in [(50e-6_f64, 1e-5), (200e-6_f64, 2e-6)] {
        let bubble = KellerMiksisBubble::try_new(
            liquid,
            // Quasistatic gas inventory at the fixed interface temperature.
            // The dynamic exponent below comes from the gas heat equation.
            BubbleGas::equilibrium(liquid, radius, reference_pressure, 1.0)?,
            KellerMiksisLimits {
                bubble: BubbleLimits {
                    minimum_radius_m: 0.5 * radius,
                    maximum_radius_m: 2.0 * radius,
                    maximum_wall_mach_number: 0.01,
                    maximum_vapor_density_ratio: 0.01,
                },
                maximum_pressure_ratio: 0.01,
            },
        )?;
        let number = fraction / (4.0 * std::f64::consts::PI * radius.powi(3) / 3.0);
        classes.push(BubblePopulationClass::try_new(bubble, number)?);
    }
    let original = PolydisperseBubblyLiquid::new(
        water,
        BubblePopulation::try_new(liquid, classes, 2)?,
        BubblyLiquidLimits {
            maximum_void_fraction: 0.001,
            maximum_wavenumber_radius: 0.1,
            maximum_wavenumber_spacing: 1.0,
            maximum_fractional_radius: 0.01,
            maximum_neighbor_scattered_pressure_ratio: 0.1,
            maximum_multiple_scattering_parameter: 0.1,
        },
    )?;
    let controls = EquilibriumControls {
        relative_radius_tolerance: 1e-13,
        absolute_radius_tolerance_m: 1e-19,
        absolute_pressure_tolerance_pa: 1e-8,
        maximum_iterations: 100,
    };
    let omega = 1000.0;
    println!("Same isothermal nuclei at 1000 rad/s; local pressure envelope 0.01 Pa");
    println!("pressure_Pa alpha R1_um R2_um phase_speed_m_s attenuation_per_m root_residual_Pa");
    for pressure in [20_000.0, 40_000.0, reference_pressure, 180_000.0] {
        let (continued, report) =
            original.at_equilibrium_pressure(pressure, &controls, &mut || false)?;
        let response = continued.response(omega, 0.01)?;
        let k = response.wavenumber_per_m();
        let radii: Vec<_> = report
            .equilibria
            .iter()
            .filter_map(|root| root.as_ref().map(|value| value.response.radius_m * 1e6))
            .collect();
        println!(
            "{pressure:.1} {:.6e} {:.6} {:.6} {:.6} {:.6e} {:.3e}",
            continued.population().void_fraction(),
            radii[0],
            radii[1],
            omega / k.re,
            k.im,
            report.maximum_absolute_pressure_residual_pa,
        );
    }

    // Evaluate the SAME gas material at each CURRENT partial pressure. A new
    // operating pressure requires new material evaluations before preparing
    // thermal response; reference gas density must never be reused blindly.
    let air = GasSpec::dry_air_ussa1976();
    let gas_states = original
        .population()
        .classes()
        .iter()
        .map(|class| -> Result<_, Box<dyn std::error::Error>> {
            let Some(linear) = class.linearization() else {
                return Ok(None);
            };
            let pressure = class
                .bubble()
                .gas()
                .pressure_at_radius(linear.equilibrium_radius_m)?;
            Ok(Some(GasState::try_new(
                &air,
                liquid.temperature_k(),
                pressure,
            )?))
        })
        .collect::<Result<Vec<_>, _>>()?;
    let thermal = ThermalPolydisperseBubblyLiquid::new(
        &original,
        &gas_states,
        ThermalBubbleLimits {
            maximum_gas_wavenumber_radius: 0.1,
            maximum_gas_liquid_density_ratio: 0.01,
            maximum_vapor_pressure_fraction: 0.1,
        },
        &mut || false,
    )?;
    let conducting = ThermalPolydisperseBubblyLiquid::new_with_conducting_liquid(
        &original,
        &gas_states,
        ThermalBubbleLimits {
            maximum_gas_wavenumber_radius: 0.1,
            maximum_gas_liquid_density_ratio: 0.01,
            maximum_vapor_pressure_fraction: 0.1,
        },
        ConductingLiquidPopulationLimits {
            maximum_penetration_to_spacing: 0.1,
            maximum_fractional_interface_temperature: 0.01,
        },
        &mut || false,
    )?;
    println!("Dry-air and IF97-water heat diffusion at the reference liquid pressure");
    println!(
        "omega_rad_s isothermal_c fixed_wall_c conducting_liquid_c fixed_wall_attenuation conducting_liquid_attenuation gas_conduction_W_m3 liquid_conduction_W_m3 liquid_viscosity_W_m3"
    );
    for omega in [100.0, 1000.0, 10_000.0, 30_000.0] {
        let isothermal = original.response(omega, 0.01)?;
        let fixed_wall = thermal.response(omega, 0.01)?;
        let response = conducting.response(omega, 0.01)?;
        let k = response.acoustic().wavenumber_per_m();
        let heat = response.mean_absorption_w_m3(0.01)?;
        println!(
            "{omega:.1} {:.6} {:.6} {:.6} {:.6e} {:.6e} {:.6e} {:.6e} {:.6e}",
            omega / isothermal.wavenumber_per_m().re,
            omega / fixed_wall.acoustic().wavenumber_per_m().re,
            omega / k.re,
            fixed_wall.acoustic().wavenumber_per_m().im,
            k.im,
            heat.gas_thermal_w_m3,
            heat.liquid_thermal_w_m3,
            heat.liquid_viscous_w_m3,
        );
    }
    Ok(())
}

// Each original member supplies its own joint nonlinear source history.
// The native acoustic renderer evaluates these at DIFFERENT retarded times
// before summing signed pressures on a single receiver clock.
fn cluster_noise(water: &SaturatedWaterCavitation) -> Result<(), Box<dyn std::error::Error>> {
    let radius = 50e-6;
    let ambient = 4000.0;
    let liquid = water.state();
    let bubble = RayleighPlessetBubble::try_new(
        liquid,
        BubbleGas::equilibrium(liquid, radius, ambient, 1.4)?,
        BubbleLimits {
            minimum_radius_m: 0.5 * radius,
            maximum_radius_m: 2.0 * radius,
            maximum_wall_mach_number: 0.01,
            maximum_vapor_density_ratio: 0.01,
        },
    )?;
    let cluster = BubbleCluster::try_new(
        vec![
            ClusterBubble {
                bubble,
                center_m: [0.0; 3],
            },
            ClusterBubble {
                bubble,
                center_m: [1e-3, 0.0, 0.0],
            },
        ],
        vec![
            BubbleState::try_new(0.0, 0.97 * radius, 0.0)?,
            BubbleState::try_new(0.0, 1.02 * radius, 0.0)?,
        ],
        0.0,
        ambient,
        BubbleClusterLimits {
            maximum_bubbles: 2,
            maximum_radius_to_separation: 0.1,
            maximum_inertia_row_sum: 0.5,
            maximum_liquid_pressure_ratio: 0.01,
            maximum_source_angular_frequency_rad_s: std::f64::consts::TAU * 20_000.0,
            maximum_retardation_phase_rad: 0.1,
        },
        ambient,
    )?;
    let mut owner = MaterialBubbleCluster::new(water, cluster)?;
    let receiver = [-0.03, 0.0, 0.0];
    let observer = owner.monopole_observer(
        receiver,
        MonopoleLimits {
            maximum_ka: 0.1,
            minimum_distance_radii: 10.0,
            maximum_pressure_ratio: 0.01,
        },
    )?;
    let numerics = BubbleClusterNumerics {
        radial: BubbleNumerics {
            relative_tolerance: 1e-9,
            radius_absolute_tolerance_m: 1e-15,
            velocity_absolute_tolerance_m_s: 1e-10,
            initial_step_s: 1e-8,
            minimum_step_s: 1e-12,
            maximum_step_s: 1e-7,
            maximum_attempts: 10_000,
        },
        absolute_energy_tolerance_j: 1e-16,
    };
    let mut histories: Vec<_> = owner
        .compact_source_samples(ambient)?
        .into_iter()
        .map(|sample| vec![sample])
        .collect();
    for index in 1..=80 {
        let step = owner.advance_observed(
            f64::from(index) * 1e-6,
            &numerics,
            &|_| ambient,
            &observer,
            &mut || false,
        )?;
        for (history, sample) in histories.iter_mut().zip(step.sources) {
            history.push(sample);
        }
    }
    let tracks: Vec<_> = histories
        .iter()
        .zip(owner.source_motions()?)
        .map(|(samples, motion)| CompactSourceHistory {
            motion,
            samples,
            source_bandwidth_hz: owner.source_bandwidth_hz(),
        })
        .collect();
    let microphone = MovingCompactObserver::new(
        liquid.acoustic_state()?,
        receiver,
        MovingCompactLimits {
            maximum_mach: 0.0,
            minimum_distance_m: 0.01,
        },
    )?;
    let mut renderer = MovingCompactRenderer::new(
        microphone,
        &tracks,
        MovingPressureConfig {
            first_receiver_time_s: 0.031 / liquid.liquid_sound_speed_m_s() + 1e-6,
            sample_rate_hz: 500_000,
            samples: 32,
            max_block: 32,
            receiver_bandwidth_hz: owner.source_bandwidth_hz(),
            maximum_time_residual_s: 1e-15,
        },
        MovingPressureBudget {
            maximum_sources: 2,
            maximum_source_samples: 162,
            maximum_output_samples: 32,
            maximum_payload_bytes: 100_000,
            maximum_source_receiver_pairs: 64,
        },
        &mut || false,
    )?;
    let mut pressure = [0.0; 32];
    renderer.block(&mut pressure)?;
    let peak = pressure.into_iter().map(f64::abs).fold(0.0, f64::max);
    println!("Two restrained interacting bubbles; original joint volume sources");
    println!(
        "receiver at x=-0.03 m: 32 coherent pressure samples at 500 kHz; peak |p'|={peak:.6e} Pa",
    );
    println!("Per-bubble retardation; frozen liquid; passive compact observation");
    Ok(())
}

// Select an explicitly isothermal gas response for the static pressure
// continuation. The dynamic growth below uses that SAME prescribed gas law;
// this example does not solve heat diffusion or gas exchange with the liquid.
fn nucleus_inception(water: &SaturatedWaterCavitation) -> Result<(), Box<dyn std::error::Error>> {
    let liquid = water.state();
    let reference_radius = 50e-6;
    let model = RayleighPlessetBubble::try_new(
        liquid,
        BubbleGas::equilibrium(liquid, reference_radius, 4000.0, 1.0)?,
        BubbleLimits {
            minimum_radius_m: 0.1 * reference_radius,
            maximum_radius_m: 10.0 * reference_radius,
            maximum_wall_mach_number: 0.01,
            maximum_vapor_density_ratio: 0.01,
        },
    )?;
    let threshold = blake_threshold(&model)?;
    let controls = EquilibriumControls {
        relative_radius_tolerance: 1e-10,
        absolute_radius_tolerance_m: 1e-15,
        absolute_pressure_tolerance_pa: 1e-8,
        maximum_iterations: 128,
    };
    // Between the fold and saturation there are two static equilibria for
    // this same nucleus. Surface tension and retained gas set their stability.
    let pressure =
        0.5 * (threshold.critical_absolute_pressure_pa + liquid.saturation_pressure_pa());
    println!("Static nucleus inception with one retained isothermal gas inventory");
    println!(
        "Blake threshold: p={:.6} Pa, radius={:.6} um",
        threshold.critical_absolute_pressure_pa,
        threshold.critical_radius_m * 1e6,
    );
    for branch in [EquilibriumBranch::Stable, EquilibriumBranch::Unstable] {
        let equilibrium = equilibrium_at_pressure(&model, pressure, branch, &controls)?;
        println!(
            "{branch:?} equilibrium below saturation: p={pressure:.6} Pa, radius={:.6} um, restoring pressure={:.6} Pa",
            equilibrium.response.radius_m * 1e6,
            equilibrium.response.restoring_pressure_pa,
        );
    }

    let growth_pressure = threshold.critical_absolute_pressure_pa - 100.0;
    match equilibrium_at_pressure(
        &model,
        growth_pressure,
        EquilibriumBranch::Stable,
        &controls,
    ) {
        Err(InceptionError::EquilibriumOutsideDomain { .. }) => {}
        Err(error) => return Err(error.into()),
        Ok(_) => {
            return Err("below-threshold nucleus unexpectedly has a static equilibrium".into());
        }
    }
    let mut state = BubbleState::try_new(0.0, threshold.critical_radius_m, 0.0)?;
    model.advance(
        &mut state,
        20e-6,
        &BubbleNumerics {
            relative_tolerance: 1e-9,
            radius_absolute_tolerance_m: 1e-14,
            velocity_absolute_tolerance_m_s: 1e-9,
            initial_step_s: 1e-8,
            minimum_step_s: 1e-12,
            maximum_step_s: 1e-7,
            maximum_attempts: 10_000,
        },
        &|_| growth_pressure,
    )?;
    println!(
        "After 20 us at p={growth_pressure:.6} Pa: radius={:.6} um, outward wall speed={:.6e} m/s",
        state.radius_m() * 1e6,
        state.wall_velocity_m_s(),
    );
    Ok(())
}
