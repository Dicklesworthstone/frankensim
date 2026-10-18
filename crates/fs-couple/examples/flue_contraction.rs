//! Fixed-window separation -> original jet/bore feedback -> exterior pressure.
//! Usage: flue_contraction NEW_OUTPUT.wav AREA_RATIO_OR_none AMBIENT_PA FULL_SCALE_PA [WINDWAY_LENGTH_M]
//! The jet card and coefficient are reduced-model data, not a calibrated flute.
use fs_aeroac::jetcard::mint_tonal_interim_card;
use fs_couple::bernoulli_aperture::contraction::ContractionSpec;
use fs_couple::bernoulli_aperture::fractional::FractionalTubePolicy;
use fs_couple::bernoulli_aperture::network::{NetworkNode,TubeNetworkSpec,TubeSection};
use fs_couple::bernoulli_aperture::radiation::BaffledRadiationLoad;
use fs_couple::flue_loop::{FlueGeometry,FlueVoice};
use fs_couple::flue_loop::network::FlueNetworkSpec;
use fs_couple::flue_loop::network::performance::{FlueOutletConfig,FlueOutletPerformance,FlueRadiationInput};
use fs_couple::pcm_wav::{encode_pcm16_wav,baffled::CircularOutletReceiver,observation::PressureRenderer};
use fs_exec::CancelGate;
use fs_material::gas::GasState;
use fs_scenario::gesture::GestureSchedule;
use std::io::Write;

fn run()->Result<(),String> {
    let args:Vec<_>=std::env::args().skip(1).collect();
    let (destination,ratio,ambient,scale,windway)=match args.as_slice() {
        [d,r,a,s] => (d,r,a,s,None),
        [d,r,a,s,l] => (d,r,a,s,Some(l)),
        _ => return Err("usage: flue_contraction NEW_OUTPUT.wav AREA_RATIO_OR_none AMBIENT_PA FULL_SCALE_PA [WINDWAY_LENGTH_M]".into()),
    };
    let parse=|x:&str|x.parse::<f64>().map_err(|_|format!("invalid SI number: {x}"));
    let ratio=if ratio=="none" {None} else {Some(parse(ratio)?)};
    let ambient=parse(ambient)?;let scale=parse(scale)?;
    let windway=windway.map(|length|parse(length)).transpose()?;
    if !scale.is_finite() || scale<=0.0 {return Err("pressure full-scale must be positive finite Pa".into());}
    let gate=CancelGate::new_clock_free();
    let mechanics=96000;let rate=48000;let samples=4800;
    let gas=GasState::try_new_moist_air(293.15,ambient,0.0).map_err(|e|format!("{e:?}"))?;
    let geometry=FlueGeometry {flue_width_m:0.01,flue_height_m:0.001,cut_up_m:0.004,
        labium_offset_m:0.0,jet_angle_rad:0.0,
        bore_area_m2:core::f64::consts::PI*0.006*0.006,mouth_loss_rel:0.02};
    let radiation=BaffledRadiationLoad::new(0.006,gas.density,gas.sound_speed,
        1.0/f64::from(mechanics),2000.0,&gate).map_err(|e|e.to_string())?;
    let source=FlueNetworkSpec {
        graph:TubeNetworkSpec {nodes:vec![NetworkNode::Inlet,radiation.termination()],
            sections:vec![TubeSection {nodes:[0,1],length_m:0.3,radius_m:0.006,max_length_error_m:1e-10}],
            sound_speed_m_s:gas.sound_speed,max_wave_memory_bytes:1<<20},
        gases:vec![gas],viscothermal:vec![],minimum_jet_speed_m_s:1.0,
        maximum_jet_mach:0.3,max_jet_history_samples:8192};
    let policy=FractionalTubePolicy {maximum_frequency_hz:2000.0,maximum_summed_phase_error_rad:0.002};
    let mut voice=FlueVoice::from_fractional_network(Some(&mint_tonal_interim_card().map_err(|e|format!("{e:?}"))?),
        90.0,1e-6,geometry,source,None,policy,mechanics,3,0,&gate).map_err(|e|e.to_string())?;
    if let Some(length)=windway {
        voice=voice.with_viscous_windway(length).map_err(|e|e.to_string())?;
    }
    if let Some(area_ratio)=ratio {
        voice=voice.with_mouth_contraction(ContractionSpec {area_ratio,maximum_contracted_mach:0.2,
            provenance:"authored symmetric fixed-mouth separation experiment; not measured flute edges".into()})
            .map_err(|e|e.to_string())?;
    }
    let schedule=GestureSchedule::from_canonical_bytes(include_bytes!("flue_contraction.gesture"))
        .map_err(|e|e.to_string())?;
    let config=FlueOutletConfig {sample_rate_hz:rate,samples,max_block:240,max_controls:256,
        max_compile_work:100000,outlet_node:1,maximum_abs_pressure_pa:1000.0,
        receiver:CircularOutletReceiver {position_m:[0.0,0.0,0.25],radial_rings:4,
            angular_points:16,maximum_frequency_hz:2000.0}};
    // Source runs at 96 kHz. Every solved outlet flow enters the existing
    // differentiating/anti-alias receiver; pressure is emitted at 48 kHz.
    let mut performance=FlueOutletPerformance::with_multirate_radiation(voice,schedule,config,
        FlueRadiationInput {exterior:gas,additional_outlets:vec![]},mechanics)
        .map_err(|e|e.to_string())?;
    let mut pressure=vec![0.0;samples as usize];
    for block in pressure.chunks_mut(240) {performance.block(block).map_err(|e|e.to_string())?;}
    let voice=performance.voice();
    let loss=voice.diagnostics().iter().map(|d|d.separation_loss_j).sum::<f64>();
    let defect=voice.diagnostics().iter().map(|d|d.ledger_defect_j.abs()).fold(0.0_f64,f64::max);
    let (wav,clips)=encode_pcm16_wav(&pressure,rate,scale).map_err(|e|e.to_string())?;
    // Do not create or overwrite output before all source/receiver/PCM checks.
    let mut file=std::fs::OpenOptions::new().write(true).create_new(true).open(destination)
        .map_err(|e|e.to_string())?;
    file.write_all(&wav).and_then(|()|file.flush()).map_err(|e|e.to_string())?;
    eprintln!("rho_kg_m3={:.9}; c_m_s={:.9}; area_ratio={ratio:?}; separation_j={loss:.9e}; max_mouth_balance_j={defect:.9e}; clips={clips}; full_scale_pa={scale}",gas.density,gas.sound_speed);
    eprintln!("windway_length_m={windway:?}; gas_viscosity_pa_s={:.9e}",gas.dynamic_viscosity);
    eprintln!("retained_network_energy_j={:.9e}; source_samples={}; output_samples={samples}; no normalization",
        voice.network().expect("physical source retained").stored_energy_j(),voice.samples_rendered());
    Ok(())
}
fn main() {
    if let Err(error)=run() {eprintln!("flue contraction refused: {error}");std::process::exit(1);}
}
