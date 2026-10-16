//! Synthetic, explicit moving-aperture/tube experiment; NOT measured reed data.
//! cargo run -p fs-couple --example aperture_tube > aperture-tube.csv
//! CSV contains interior/terminal pressures, not a far-field microphone signal.

use fs_couple::bernoulli_aperture::BernoulliAperture;
use fs_couple::bernoulli_aperture::dynamic::{
    ApertureState, ApertureTerminal, DynamicAperture, DynamicApertureSpec,
};
use fs_couple::bernoulli_aperture::tube::{ApertureTube, TubeDrive, TubeFrame, UniformTubeSpec};
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;
use std::io::Write;

fn run() -> Result<(), String> {
    let dt = 1e-5;
    let tube = UniformTubeSpec {
        length_m: 0.25, radius_m: 0.007, sound_speed_m_s: 343.0,
        terminal_reflection: -0.8, max_length_error_m: 0.0005,
        max_wave_memory_bytes: 1 << 20,
    };
    let mechanics = DynamicApertureSpec {
        aperture: BernoulliAperture {
            rest_opening_m: 4e-4, width_m: 0.013, closing_pressure_pa: 6000.0,
        },
        mass_kg: 1e-5, stiffness_n_m: 500.0, damping_ratio: 0.35,
        density_kg_m3: 1.2,
        impedance_pa_s_m3: tube.characteristic_impedance(1.2).map_err(|e| format!("{e:?}"))?,
        time_step_s: dt, max_steps: 4096,
    };
    let contact = Obstacle::new(vec![-1.0], 1, 1, vec![0.0], vec![1.0], 1e8, 2.0,
        "synthetic example, not identified constitutive data".into())
        .and_then(|o| o.with_internal_loss(5.0)).map_err(|e| e.to_string())?;
    let aperture = DynamicAperture::new(mechanics, ApertureState {
        opening_m: mechanics.aperture.rest_opening_m, opening_velocity_m_s: 0.0,
    }, contact).map_err(|e| format!("{e:?}"))?;
    let mut model = ApertureTube::new(aperture, tube).map_err(|e| format!("{e:?}"))?;
    eprintln!("synthetic_model; requested_length_m={}; represented_length_m={}; one_way_samples={}",
        tube.length_m, model.represented_length_m(), model.one_way_samples());
    let stdout = std::io::stdout();
    let mut out = std::io::BufWriter::new(stdout.lock());
    writeln!(out, "step_end_time_s,midpoint_bore_pa,terminal_pa,opening_m,total_energy_j,balance_residual_j")
        .map_err(|e| e.to_string())?;
    let gate = CancelGate::new_clock_free();
    let mut inputs = [TubeDrive::default(); 128];
    let mut frames = [TubeFrame::default(); 128];
    for block in 0..32 {
        for (i, drive) in inputs.iter_mut().enumerate() {
            *drive = TubeDrive {
                upstream_pressure_pa: if block * 128 + i < 1024 { 1200.0 } else { 0.0 },
                body_flow_m3_s: 0.0,
            };
        }
        let progress = model.advance_block(&inputs, &mut frames, &gate).map_err(|e| format!("{e:?}"))?;
        for frame in &frames[..progress.completed] {
            writeln!(out, "{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
                frame.aperture.time_s, frame.aperture.bore_pressure_pa,
                frame.waveguide.terminal_pressure_pa, frame.aperture.state.opening_m,
                frame.stored_energy_j, frame.balance_residual_j()).map_err(|e| e.to_string())?;
        }
        if progress.terminal != ApertureTerminal::Complete {
            return Err(format!("stopped at accepted sample {}: {:?}",
                model.aperture().accepted_steps(), progress.terminal));
        }
    }
    out.flush().map_err(|e| e.to_string())
}

fn main() {
    if let Err(error) = run() {
        eprintln!("aperture/tube experiment failed: {error}");
        std::process::exit(1);
    }
}
