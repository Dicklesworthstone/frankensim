//! Synthetic locally reacting wall experiment; NOT identified material data.
//! cargo run -p fs-couple --example aperture_wall > aperture-wall.csv
//! Pressures are inside the tube, not a predicted far-field microphone signal.

use fs_couple::bernoulli_aperture::BernoulliAperture;
use fs_couple::bernoulli_aperture::dynamic::{ApertureState, DynamicAperture, DynamicApertureSpec};
use fs_couple::bernoulli_aperture::network::ApertureNetwork;
use fs_couple::bernoulli_aperture::tube::{TubeDrive, UniformTubeSpec};
use fs_couple::bernoulli_aperture::wall::{WallPin, lined_tube};
use fs_dcontact::Obstacle;
use std::io::Write;

fn run() -> Result<(), String> {
    let dt = 1e-5;
    let tube = UniformTubeSpec {
        length_m: 64.0 * 343.0 * dt, radius_m: 0.007, sound_speed_m_s: 343.0,
        terminal_reflection: -0.8, max_length_error_m: 1e-12, max_wave_memory_bytes: 1 << 20,
    };
    let lined = lined_tube(tube, WallPin {
        surface_density: 0.2, stiffness_per_area: 2e7, resistance: 300.0,
    }, 4, 1 << 20).map_err(|e| format!("{e:?}"))?;
    let mechanics = DynamicApertureSpec {
        aperture: BernoulliAperture { rest_opening_m: 4e-4, width_m: 0.013, closing_pressure_pa: 6000.0 },
        mass_kg: 1e-5, stiffness_n_m: 500.0, damping_ratio: 0.35,
        density_kg_m3: 1.2, impedance_pa_s_m3: lined.network.inlet_impedance(1.2).map_err(|e| format!("{e:?}"))?,
        time_step_s: dt, max_steps: 4096,
    };
    let contact = Obstacle::new(vec![-1.0], 1, 1, vec![0.0], vec![1.0], 1e8, 2.0,
        "synthetic wall experiment contact; not identified data".into())
        .and_then(|o| o.with_internal_loss(5.0)).map_err(|e| e.to_string())?;
    let aperture = DynamicAperture::new(mechanics, ApertureState {
        opening_m: mechanics.aperture.rest_opening_m, opening_velocity_m_s: 0.0,
    }, contact).map_err(|e| format!("{e:?}"))?;
    let mut model = ApertureNetwork::new(aperture, lined.network).map_err(|e| format!("{e:?}"))?;
    eprintln!("synthetic_model; wall_patches={}; requested_length_m={}; represented_length_m={}",
        lined.patches.len(), tube.length_m,
        model.represented_sections().iter().map(|s| s.represented_length_m).sum::<f64>());
    let stdout = std::io::stdout();
    let mut out = std::io::BufWriter::new(stdout.lock());
    writeln!(out, "step_end_time_s,midpoint_inlet_pa,valve_opening_m,first_wall_displacement_m,first_wall_velocity_m_s,wall_energy_j,wave_energy_j,total_energy_j,wall_loss_j,balance_residual_j")
        .map_err(|e| e.to_string())?;
    for sample in 0..4096 {
        let frame = model.step(TubeDrive {
            upstream_pressure_pa: if sample < 1024 { 1200.0 } else { 0.0 }, body_flow_m3_s: 0.0,
        }).map_err(|e| format!("{e:?}"))?;
        // Read immediately after this accepted sample, not after an entire block:
        // the exported wall coordinates must belong to the exported frame.
        let wall = lined.patches[0].observe(&model, 1).map_err(|e| format!("{e:?}"))?;
        writeln!(out, "{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
            frame.aperture.time_s, frame.aperture.bore_pressure_pa, frame.aperture.state.opening_m,
            wall.displacement_m, wall.velocity_m_s, frame.network.load_stored_energy_j,
            frame.network.wave_stored_energy_j, frame.stored_energy_j,
            frame.network.interior_loss_j, frame.balance_residual_j()).map_err(|e| e.to_string())?;
    }
    out.flush().map_err(|e| e.to_string())
}
fn main() {
    if let Err(error) = run() {
        eprintln!("compliant-wall experiment failed: {error}");
        std::process::exit(1);
    }
}
