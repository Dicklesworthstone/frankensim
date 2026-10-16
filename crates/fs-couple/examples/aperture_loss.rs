//! Synthetic frequency-dependent cavity loss, not measured material data.
//! cargo run -p fs-couple --example aperture_loss > aperture-loss.csv
//! Pressure is internal, not a far-field microphone or a calibrated waveform.
use fs_couple::bernoulli_aperture::BernoulliAperture;
use fs_couple::bernoulli_aperture::cavity::HelmholtzLoadSpec;
use fs_couple::bernoulli_aperture::dynamic::{ApertureState, DynamicAperture, DynamicApertureSpec};
use fs_couple::bernoulli_aperture::loss::fit_boundary_loss;
use fs_couple::bernoulli_aperture::network::{ApertureNetwork, NetworkNode, TubeNetworkSpec, TubeSection};
use fs_couple::bernoulli_aperture::tube::TubeDrive;
use fs_dcontact::Obstacle;
use std::io::Write;

fn run() -> Result<(), String> {
    let cavity = HelmholtzLoadSpec { volume_m3: 1e-4, neck_radius_m: 0.003,
        effective_neck_length_m: 0.02, resistance_pa_s_m3: 2e5 };
    let base = cavity.impedance(1.2, 343.0).map_err(|e| format!("{e:?}"))?;
    // Explicit synthetic EXCESS resistance observations, separate from the
    // constant 2e5 base resistance. Check points were not used for fitting.
    let target = |w: f64| 2e5 * w * w / (w * w + 500.0_f64.powi(2))
        + 3e5 * w * w / (w * w + 6000.0_f64.powi(2));
    let training = [500.0, 6000.0].map(|w| (w, target(w)));
    let checks = [1000.0, 2000.0, 4000.0].map(|w| (w, target(w)));
    let fit = fit_boundary_loss(base, &training, &checks, 1e-7, 1e-10)
        .map_err(|e| format!("{e:?}"))?;
    eprintln!("synthetic_excess_loss; checked_error_pa_s_m3={:.17e}; poles={:?}; no_continuous_band_or_imaginary_part_validation",
        fit.max_checked_error_pa_s_m3, fit.load.terms());
    let sections = [(0, 1, 0.12, 0.007), (1, 2, 0.18, 0.009), (1, 3, 0.02, 0.003)]
        .map(|(a, b, length, radius)| TubeSection { nodes: [a, b], length_m: length,
            radius_m: radius, max_length_error_m: 0.0018 }).to_vec();
    let network = TubeNetworkSpec { nodes: vec![NetworkNode::Inlet, NetworkNode::Junction,
        NetworkNode::Termination { reflection: -0.8 }, fit.termination()], sections,
        sound_speed_m_s: 343.0, max_wave_memory_bytes: 1 << 20 };
    let mechanics = DynamicApertureSpec {
        aperture: BernoulliAperture { rest_opening_m: 4e-4, width_m: 0.013, closing_pressure_pa: 6000.0 },
        mass_kg: 1e-5, stiffness_n_m: 500.0, damping_ratio: 0.35, density_kg_m3: 1.2,
        impedance_pa_s_m3: network.inlet_impedance(1.2).map_err(|e| format!("{e:?}"))?,
        time_step_s: 1e-5, max_steps: 4096,
    };
    let lay = Obstacle::new(vec![-1.0], 1, 1, vec![0.0], vec![1.0], 1e8, 2.0,
        "synthetic example contact law".into()).and_then(|l| l.with_internal_loss(5.0))
        .map_err(|e| e.to_string())?;
    let aperture = DynamicAperture::new(mechanics, ApertureState {
        opening_m: 4e-4, opening_velocity_m_s: 0.0 }, lay).map_err(|e| format!("{e:?}"))?;
    let mut model = ApertureNetwork::new(aperture, network).map_err(|e| format!("{e:?}"))?;
    let stdout = std::io::stdout();
    let mut out = std::io::BufWriter::new(stdout.lock());
    writeln!(out, "time_s,inlet_pa,cavity_port_pa,opening_m,load_energy_j,total_energy_j,dissipation_j,balance_residual_j")
        .map_err(|e| e.to_string())?;
    for n in 0..4096 {
        let f = model.step(TubeDrive { upstream_pressure_pa: if n < 1024 { 1200.0 } else { 0.0 },
            body_flow_m3_s: 0.0 }).map_err(|e| format!("{e:?}"))?;
        let p = model.node_frame(3).expect("declared cavity terminal").pressure_pa;
        writeln!(out, "{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
            f.aperture.time_s, f.aperture.bore_pressure_pa, p, f.aperture.state.opening_m,
            f.network.load_stored_energy_j, f.stored_energy_j, f.dissipated_energy_j, f.balance_residual_j())
            .map_err(|e| e.to_string())?;
    }
    out.flush().map_err(|e| e.to_string())
}
fn main() {
    if let Err(e) = run() { eprintln!("aperture loss example failed: {e}"); std::process::exit(1); }
}
