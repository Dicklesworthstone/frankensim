//! Helium-II pressure/entropy waves driving an actual native finite-mass piston.
//!
//! Run with `cargo run -p fs-couple --features superfluid --example superfluid_piston`.
//! An optional first argument changes the cell count (default 64).
//! Helium properties use the sourced 1.8 K saturated-vapor-pressure reference.
//! Geometry, excitation, piston dimensions and mutual-friction coefficient are
//! authored demonstration inputs, not a calibrated cryogenic apparatus.
//! Reference pressure is balanced outside the piston. Motion is a small
//! displacement on the original fixed-reference channel; exported heat is
//! accounted for but this example does not feed it back into that reference.

use fs_couple::superfluid::piston::{
    BalancedPistonBoundary, PistonControls, advance_balanced_piston,
};
use fs_duct::superfluid::{ChannelEnd, LinearizationLimits, MutualFriction, TwoFluidChannel};
use fs_material::superfluid::helium::HeliumSvpReference;
use fs_mbd::{MassProperties, Pose, RigidBodyState, Vec3};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cells = std::env::args().nth(1).map_or(Ok(64usize), |s| s.parse())?;
    if !(8..=4096).contains(&cells) {
        return Err("cell count must be between 8 and 4096".into());
    }
    let reference = HeliumSvpReference::at_1_8_k()?;
    let material = *reference.state();
    let length_m = 0.1;
    let area_m2 = 1e-4;
    let mut channel = TwoFluidChannel::try_new(
        material,
        length_m,
        area_m2,
        cells,
        LinearizationLimits {
            max_pressure_perturbation_pa: 10.0,
            max_temperature_perturbation_k: 1e-3,
            max_component_speed_m_s: 0.1,
        },
    )?;
    let pressure: Vec<_> = (0..cells)
        .map(|i| fs_math::det::cos(core::f64::consts::PI * (i as f64 + 0.5) / cells as f64))
        .collect();
    let temperature: Vec<_> = (0..cells)
        .map(|i| {
            2e-5 * fs_math::det::cos(2.0 * core::f64::consts::PI * (i as f64 + 0.5) / cells as f64)
        })
        .collect();
    let still = vec![0.0; cells + 1];
    channel.initialize(&pressure, &temperature, &still, &still, || false)?;

    // A centered 2 g disk, radius 5 mm and thickness 1 mm. No spin is imposed.
    let piston_mass_kg = 0.002;
    let radius_m: f64 = 0.005;
    let thickness_m: f64 = 0.001;
    let axial_inertia = 0.5 * piston_mass_kg * radius_m * radius_m;
    let transverse_inertia =
        piston_mass_kg * (3.0 * radius_m * radius_m + thickness_m * thickness_m) / 12.0;
    let properties = MassProperties::new(
        piston_mass_kg,
        Vec3::ZERO,
        Vec3::new(axial_inertia, transverse_inertia, transverse_inertia),
    )?;
    let mut body = RigidBodyState::new(Pose::identity(), Vec3::ZERO, Vec3::ZERO)?;
    let boundary = BalancedPistonBoundary {
        end: ChannelEnd::Right,
        positive_axis_world: Vec3::new(1.0, 0.0, 0.0),
        equilibrium_position_world: Vec3::ZERO,
        max_displacement_m: 0.001,
    };
    let controls = PistonControls {
        momentum_tolerance_n_s: 1e-15,
        energy_tolerance_j: 1e-19,
        velocity_tolerance_m_s: 1e-12,
    };
    let duration_s = 0.005;
    let steps = (duration_s / (0.8 * channel.maximum_step_s())).ceil() as usize;
    let dt_s = duration_s / steps as f64;
    let initial_energy = channel.wave_energy_j()?;
    let initial_momentum = channel.resolved_momentum_kg_m_s()?;
    let half_mass = 0.5 * material.density_kg_m3() * area_m2 * length_m / cells as f64;
    let mut support_impulse = 0.0;
    let mut heat = 0.0;
    for _ in 0..steps {
        let report = advance_balanced_piston(
            &mut channel,
            &mut body,
            properties,
            boundary,
            dt_s,
            MutualFriction::Linear {
                gamma_kg_m3_s: 1000.0,
            },
            controls,
            || false,
        )?;
        support_impulse += report.stationary_support_impulse_x_n_s;
        heat += report.channel.mutual_friction_heat_j;
    }
    let velocity = body.center_of_mass_velocity_world(properties)?.x;
    let energy_defect = channel.wave_energy_j()?
        + body.kinetic_energy(properties)?
        + 0.5 * half_mass * velocity * velocity
        + heat
        - initial_energy;
    let momentum_defect = channel.resolved_momentum_kg_m_s()?
        + body.linear_momentum_world().x
        + half_mass * velocity
        + support_impulse
        - initial_momentum;
    if energy_defect.abs() > 1e-8 * initial_energy || momentum_defect.abs() > 1e-13 {
        return Err("cumulative physical momentum or energy balance failed".into());
    }
    println!(
        "helium_reference_temperature_k={:.9}",
        material.temperature_k()
    );
    println!("cells={cells} steps={steps}");
    println!("simulated_time_s={:.12e}", channel.time_s());
    println!("fast_sound_m_s={:.9}", material.sound_modes().fast_m_s);
    println!("slow_sound_m_s={:.9}", material.sound_modes().slow_m_s);
    println!("piston_mass_kg={piston_mass_kg:.12e}");
    println!("boundary_fluid_half_mass_kg={half_mass:.12e}");
    println!(
        "piston_displacement_x_m={:.12e}",
        body.pose().position_world().x
    );
    println!("piston_velocity_x_m_s={velocity:.12e}");
    println!("exported_mutual_friction_heat_j={heat:.12e}");
    println!(
        "relative_combined_energy_defect={:.12e}",
        energy_defect / initial_energy
    );
    println!("combined_momentum_defect_n_s={momentum_defect:.12e}");
    Ok(())
}
