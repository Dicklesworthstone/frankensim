//! A measured 1.8 K helium-II reference drives a closed two-fluid channel.
//! Run: cargo run -p fs-duct --features superfluid --example superfluid_second_sound
//! Optional arguments: cell count, explicitly authored gamma [kg/(m³ s)].

use fs_duct::superfluid::{
    BoundaryVelocities, LinearizationLimits, MutualFriction, TwoFluidChannel,
};
use fs_material::superfluid::helium::{HeliumSvpReference, SOURCE_DOI};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.len() > 2 {
        return Err("usage: superfluid_second_sound [cells] [gamma_kg_m3_s]".into());
    }
    let cells: usize = args.first().map_or(Ok(128), |s| s.parse())?;
    if !(8..=65536).contains(&cells) {
        return Err("this example admits 8 through 65536 cells".into());
    }
    let gamma: f64 = args.get(1).map_or(Ok(0.0), |s| s.parse())?;
    if !gamma.is_finite() || gamma < 0.0 {
        return Err("gamma must be explicitly nonnegative and finite".into());
    }
    let reference = HeliumSvpReference::at_1_8_k()?;
    let material = *reference.state();
    let sound = material.sound_modes();
    let length = 0.1;
    let mut channel = TwoFluidChannel::try_new(
        material,
        length,
        1e-4,
        cells,
        LinearizationLimits {
            max_pressure_perturbation_pa: 10.0,
            max_temperature_perturbation_k: 1e-4,
            max_component_speed_m_s: 0.01,
        },
    )?;
    let temperature: Vec<_> = (0..cells)
        .map(|i| 1e-6 * (std::f64::consts::PI * (i as f64 + 0.5) / cells as f64).cos())
        .collect();
    channel.initialize(
        &vec![0.0; cells],
        &temperature,
        &vec![0.0; cells + 1],
        &vec![0.0; cells + 1],
        || false,
    )?;
    let period = 2.0 * length / sound.slow_m_s;
    let steps = (period / (0.75 * channel.maximum_step_s())).ceil() as usize;
    let dt = period / steps as f64;
    let friction = if gamma == 0.0 {
        MutualFriction::None
    } else {
        MutualFriction::Linear {
            gamma_kg_m3_s: gamma,
        }
    };
    let initial = channel.wave_energy_j()?;
    println!(
        "# source DOI {SOURCE_DOI}; one SVP reference, T=1.8 K; p={:.6} Pa",
        reference.saturation_pressure_pa()
    );
    println!(
        "# first sound {:.9} m/s; predicted second sound {:.9} m/s; independently observed {:.2} m/s",
        sound.fast_m_s,
        sound.slow_m_s,
        reference.raw().second_sound_m_s
    );
    println!(
        "# cells={cells}; steps={steps}; dt={dt:.9e} s; authored gamma={gamma:.9e} kg/(m^3 s)"
    );
    println!(
        "# heat is exported for a thermal owner; this fixed-reference run does not update its background temperature"
    );
    println!(
        "time_s,pressure_left_pa,temperature_left_delta_k,normal_velocity_mid_m_s,superfluid_velocity_mid_m_s,counterflow_heat_flux_mid_w_m2,wave_energy_j,exported_friction_heat_j"
    );
    emit(&channel)?;
    for step in 1..=steps {
        channel.step(dt, BoundaryVelocities::CLOSED, friction, || false)?;
        if step % (steps / 8).max(1) == 0 || step == steps {
            emit(&channel)?;
        }
    }
    let residual = channel.wave_energy_j()? + channel.cumulative_mutual_friction_heat_j() - initial;
    println!(
        "# relative closed wave-plus-exported-heat defect={:.9e}",
        residual / initial
    );
    Ok(())
}

fn emit(channel: &TwoFluidChannel) -> Result<(), Box<dyn std::error::Error>> {
    let cell = channel.cell(0).ok_or("empty channel")?;
    let face = channel
        .face(channel.cells() / 2)
        .ok_or("missing middle face")?;
    println!(
        "{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e}",
        channel.time_s(),
        cell.pressure_perturbation_pa,
        cell.temperature_perturbation_k,
        face.normal_velocity_m_s,
        face.superfluid_velocity_m_s,
        face.counterflow_heat_flux_w_m2,
        channel.wave_energy_j()?,
        channel.cumulative_mutual_friction_heat_j()
    );
    Ok(())
}
