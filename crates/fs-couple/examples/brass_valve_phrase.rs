//! Manufactured two-path control demonstration, NOT a radiated instrument.
//! cargo run -p fs-couple --example brass_valve_phrase -- mouthpiece.wav
use fs_couple::brass_loop::{BrassControl, LipIsland};
use fs_couple::brass_sequence::{BrassEvent, BrassSequence, BrassSequenceConfig};
use fs_couple::driving_point::LineLadder;
use fs_couple::mm_line::{MmBankLoad, MmLineConfig, MmLoad};
use fs_couple::pcm_wav::{encode_pcm16_wav, observation::PressureRenderer};
use fs_duct::{Duct, LossModel, Segment, TabulatedLoad, TabulatedMouth};
use fs_material::gas::{GasSpec, GasState};
use fs_math::c64::C64;
use std::io::Write;

fn run() -> Result<(), String> {
    let args: Vec<_> = std::env::args_os().collect();
    if args.len() != 2 { return Err("usage: brass_valve_phrase NEW_OUTPUT.wav".into()); }
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0)
        .map_err(|e| format!("gas: {e:?}"))?;
    let rate = 48000;
    let dt = 1.0 / f64::from(rate);
    let radius = 0.002;
    let zc = gas.density * gas.sound_speed / (core::f64::consts::PI * radius * radius);
    let mut ducts = Vec::new();
    let mut tables = Vec::new();
    for (transit, resistance) in [(8.0, 0.7), (12.0, 1.3)] {
        let length = transit * gas.sound_speed * dt;
        ducts.push(Duct { segments: vec![Segment::Cylinder { radius, length }] });
        tables.push(TabulatedLoad::try_new(
            vec![core::f64::consts::TAU / (LineLadder::default().max_taps as f64 * dt),
                core::f64::consts::PI / dt],
            vec![C64::new(resistance * zc, 0.0); 2],
            "manufactured resistive mouth, not measured radiation", 0.0,
        ).map_err(|e| e.to_string())?.with_mouth(TabulatedMouth {
            axial_position_m: length, radius_m: radius, gas,
            dc_impedance: resistance * zc,
        }).map_err(|e| e.to_string())?);
    }
    let loads: Vec<_> = tables.iter().map(MmLoad::Tabulated).collect();
    let line = MmLineConfig { sample_rate_hz: rate, n_modes: 2,
        extra_slices: 1, loss: LossModel::Lossless };
    let lip = LipIsland { mass_kg: 1.8e-4, stiffness_n_m: 350.0, damping_n_s_m: 0.05,
        width_m: 0.012, rest_gap_m: 5e-4, face_area_m2: 2e-5,
        provenance: "manufactured control example; no specimen qualification".into() };
    let events = vec![
        BrassEvent { sample: 1200, control: BrassControl::SetBlowingPressure(1500.0) },
        BrassEvent { sample: 4801, control: BrassControl::SetValve { combo: 1, fade_samples: 96 } },
        BrassEvent { sample: 9607, control: BrassControl::SetLipTension(1.2) },
        BrassEvent { sample: 14401, control: BrassControl::SetValve { combo: 0, fade_samples: 96 } },
        BrassEvent { sample: 19200, control: BrassControl::SetBlowingPressure(0.0) },
    ];
    let config = BrassSequenceConfig { samples: 24000, max_block_len: 257,
        initial_blowing_pressure_pa: 0.0, initial_lip_tension: 1.0, initial_combo: 0 };
    let mut sequence = BrassSequence::new(&ducts, &["short", "long"], &gas,
        MmBankLoad::PerCombo(&loads), &line, lip, 1.5e-6, config, events)
        .map_err(|e| format!("sequence: {e:?}"))?;
    let mut pressure = vec![0.0; config.samples as usize];
    for chunk in pressure.chunks_mut(config.max_block_len) {
        sequence.block(chunk).map_err(|e| format!("physical callback: {e:?}"))?;
    }
    let (wav, clips) = encode_pcm16_wav(&pressure, rate, 5000.0).map_err(|e| e.to_string())?;
    let mut output = std::fs::OpenOptions::new().write(true).create_new(true)
        .open(&args[1]).map_err(|e| e.to_string())?;
    output.write_all(&wav).map_err(|e| e.to_string())?;
    println!("Mouthpiece pressure only; manufactured geometry/lip/terminal loads. samples={} rate_hz={} full_scale_pa=5000 clips={} valve_changes={}",
        pressure.len(), rate, clips, sequence.voice().lift_log().len());
    Ok(())
}
fn main() {
    if let Err(error) = run() { eprintln!("{error}"); std::process::exit(1); }
}
