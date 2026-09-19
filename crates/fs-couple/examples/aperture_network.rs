//! Synthetic branched-tube experiment, not measured instrument parameters.
//! cargo run -p fs-couple --example aperture_network > aperture-network.csv
//! Endpoint pressure is NOT a far-field microphone or a radiation prediction.

use fs_couple::bernoulli_aperture::BernoulliAperture;
use fs_couple::bernoulli_aperture::dynamic::{ApertureState, DynamicAperture, DynamicApertureSpec};
use fs_couple::bernoulli_aperture::network::{ApertureNetwork, NetworkNode, TubeNetworkSpec, TubeSection};
use fs_couple::bernoulli_aperture::tube::TubeDrive;
use fs_dcontact::Obstacle;
use std::io::Write;

fn run() -> Result<(), String> {
    let dt = 1e-5;
    let section = |a, b, length_m, radius_m| TubeSection {
        nodes: [a, b], length_m, radius_m, max_length_error_m: 0.0018,
    };
    let spec = TubeNetworkSpec {
        nodes: vec![NetworkNode::Inlet, NetworkNode::Junction,
            NetworkNode::Termination { reflection: -0.8 }, NetworkNode::Termination { reflection: 1.0 }],
        sections: vec![section(0, 1, 0.12, 0.007), section(1, 2, 0.18, 0.009), section(1, 3, 0.02, 0.003)],
        sound_speed_m_s: 343.0, max_wave_memory_bytes: 1 << 20,
    };
    let mechanics = DynamicApertureSpec {
        aperture: BernoulliAperture { rest_opening_m: 4e-4, width_m: 0.013, closing_pressure_pa: 6000.0 },
        mass_kg: 1e-5, stiffness_n_m: 500.0, damping_ratio: 0.35,
        density_kg_m3: 1.2, impedance_pa_s_m3: spec.inlet_impedance(1.2).map_err(|e| format!("{e:?}"))?,
        time_step_s: dt, max_steps: 4096,
    };
    let lay = Obstacle::new(vec![-1.0], 1, 1, vec![0.0], vec![1.0], 1e8, 2.0,
        "synthetic network example, not identified material data".into())
        .and_then(|o| o.with_internal_loss(5.0)).map_err(|e| e.to_string())?;
    let valve = DynamicAperture::new(mechanics, ApertureState {
        opening_m: mechanics.aperture.rest_opening_m, opening_velocity_m_s: 0.0,
    }, lay).map_err(|e| format!("{e:?}"))?;
    let mut network = ApertureNetwork::new(valve, spec).map_err(|e| format!("{e:?}"))?;
    for (i, section) in network.represented_sections().iter().enumerate() {
        eprintln!("section={i}; requested_m={}; represented_m={}; one_way_samples={}; impedance={}",
            network.spec().sections[i].length_m, section.represented_length_m,
            section.one_way_samples, section.impedance_pa_s_m3);
    }
    let stdout = std::io::stdout();
    let mut out = std::io::BufWriter::new(stdout.lock());
    writeln!(out, "time_s,inlet_pa,junction_pa,branch_end_pa,opening_m,total_energy_j,balance_residual_j")
        .map_err(|e| e.to_string())?;
    for n in 0..4096 {
        // An explicit ideal pressure-release/closed boundary switch, not a
        // physical pad-motion model. Existing waves and valve motion persist.
        if n == 512 || n == 2048 {
            network.set_terminal_reflection(3, if n == 512 { -1.0 } else { 1.0 })
                .map_err(|e| format!("{e:?}"))?;
        }
        let frame = network.step(TubeDrive {
            upstream_pressure_pa: if n < 1024 { 1200.0 } else { 0.0 }, body_flow_m3_s: 0.0,
        }).map_err(|e| format!("{e:?}"))?;
        writeln!(out, "{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
            frame.aperture.time_s, frame.aperture.bore_pressure_pa,
            network.node_frame(1).expect("declared junction").pressure_pa,
            network.node_frame(3).expect("declared branch end").pressure_pa,
            frame.aperture.state.opening_m, frame.stored_energy_j, frame.balance_residual_j())
            .map_err(|e| e.to_string())?;
    }
    out.flush().map_err(|e| e.to_string())
}
fn main() {
    if let Err(error) = run() {
        eprintln!("aperture/network experiment failed: {error}");
        std::process::exit(1);
    }
}
