//! Synthetic moving valve / branch cavity experiment, NOT measured material data.
//! cargo run -p fs-couple --example aperture_cavity > aperture-cavity.csv
//! Pressures are internal port observations, not a far-field microphone signal.
use fs_couple::bernoulli_aperture::BernoulliAperture;
use fs_couple::bernoulli_aperture::cavity::HelmholtzLoadSpec;
use fs_couple::bernoulli_aperture::dynamic::{ApertureState,DynamicAperture,DynamicApertureSpec};
use fs_couple::bernoulli_aperture::network::{ApertureNetwork,NetworkNode,TubeNetworkSpec,TubeSection};
use fs_couple::bernoulli_aperture::tube::TubeDrive;
use fs_dcontact::Obstacle;
use std::io::Write;

fn run() -> Result<(),String> {
    let cavity=HelmholtzLoadSpec { volume_m3:1e-4,neck_radius_m:0.003,
        effective_neck_length_m:0.02,resistance_pa_s_m3:2e5 };
    let section=|a,b,length,radius| TubeSection {nodes:[a,b],length_m:length,
        radius_m:radius,max_length_error_m:0.0018};
    let spec=TubeNetworkSpec {
        nodes:vec![NetworkNode::Inlet,NetworkNode::Junction,NetworkNode::Termination {reflection:-0.8},
            cavity.termination(1.2,343.0).map_err(|e|format!("{e:?}"))?],
        // Side branch is a connection tube; the lumped cavity neck above is
        // separate physical geometry and is NOT counted a second time here.
        sections:vec![section(0,1,0.12,0.007),section(1,2,0.18,0.009),section(1,3,0.02,0.003)],
        sound_speed_m_s:343.0,max_wave_memory_bytes:1<<20,
    };
    let mechanics=DynamicApertureSpec {
        aperture:BernoulliAperture {rest_opening_m:4e-4,width_m:0.013,closing_pressure_pa:6000.0},
        mass_kg:1e-5,stiffness_n_m:500.0,damping_ratio:0.35,density_kg_m3:1.2,
        impedance_pa_s_m3:spec.inlet_impedance(1.2).map_err(|e|format!("{e:?}"))?,
        time_step_s:1e-5,max_steps:4096,
    };
    let contact=Obstacle::new(vec![-1.0],1,1,vec![0.0],vec![1.0],1e8,2.0,
        "synthetic example contact; not identified data".into())
        .and_then(|o|o.with_internal_loss(5.0)).map_err(|e|e.to_string())?;
    let aperture=DynamicAperture::new(mechanics,ApertureState {opening_m:4e-4,opening_velocity_m_s:0.0},contact)
        .map_err(|e|format!("{e:?}"))?;
    let mut net=ApertureNetwork::new(aperture,spec).map_err(|e|format!("{e:?}"))?;
    eprintln!("synthetic_lumped_cavity; continuous_unloaded_resonance_hz={:.9}; runtime_bilinear_warping=true",
        cavity.resonance_hz(1.2,343.0).map_err(|e|format!("{e:?}"))?);
    let stdout=std::io::stdout();let mut out=std::io::BufWriter::new(stdout.lock());
    writeln!(out,"time_s,inlet_pa,cavity_port_pa,opening_m,cavity_storage_j,wave_storage_j,total_energy_j,balance_j")
        .map_err(|e|e.to_string())?;
    for n in 0..4096 {
        let f=net.step(TubeDrive {upstream_pressure_pa:if n<1024 {1200.0} else {0.0},body_flow_m3_s:0.0})
            .map_err(|e|format!("sample {n}: {e:?}"))?;
        let terminal=net.node_frame(3).ok_or("missing cavity observation")?;
        writeln!(out,"{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
            f.aperture.time_s,f.aperture.bore_pressure_pa,terminal.pressure_pa,
            f.aperture.state.opening_m,f.network.load_stored_energy_j,
            f.network.wave_stored_energy_j,f.stored_energy_j,f.balance_residual_j()).map_err(|e|e.to_string())?;
    }
    out.flush().map_err(|e|e.to_string())
}
fn main() {
    if let Err(error)=run() { eprintln!("aperture/cavity experiment failed: {error}");std::process::exit(1); }
}
