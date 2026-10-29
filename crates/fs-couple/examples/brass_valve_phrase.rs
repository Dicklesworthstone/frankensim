//! Manufactured two-path control demonstration, NOT a radiated instrument.
//! cargo run -p fs-couple --example brass_valve_phrase -- mouthpiece.wav mouth.csv
use fs_couple::brass_loop::{BrassControl, LipIsland};
use fs_couple::brass_sequence::{BrassEvent, BrassSequenceConfig,
    mouth::{BrassMouthSequence, MouthFlowConfig}};
use fs_couple::driving_point::LineLadder;
use fs_couple::mm_line::{MmBankLoad, MmLineConfig, MmLoad};
use fs_couple::pcm_wav::encode_pcm16_wav;
use fs_duct::{Duct, LossModel, Segment, TabulatedLoad, TabulatedMouth};
use fs_material::gas::{GasSpec, GasState};
use fs_math::c64::C64;
use std::io::Write;

fn run() -> Result<(), String> {
    let args: Vec<_> = std::env::args_os().collect();
    if !(2..=3).contains(&args.len()) {
        return Err("usage: brass_valve_phrase NEW_OUTPUT.wav [NEW_MOUTH_FLOW.csv]".into());
    }
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
    let mut sequence = BrassMouthSequence::new(&ducts, &["short", "long"], &gas,
        MmBankLoad::PerCombo(&loads), &line, lip, 1.5e-6, config, events,
        MouthFlowConfig::default()).map_err(|e| format!("sequence: {e:?}"))?;
    let mut pressure = vec![0.0; config.samples as usize];
    let mut outgoing = pressure.clone();
    let modes = sequence.mode_count();
    let mut mouth = vec![0.0; pressure.len() * modes];
    for ((p, a), u) in pressure.chunks_mut(config.max_block_len)
        .zip(outgoing.chunks_mut(config.max_block_len))
        .zip(mouth.chunks_mut(config.max_block_len * modes))
    {
        sequence.block(p, a, u).map_err(|e| format!("physical callback: {e:?}"))?;
    }
    let (wav, clips) = encode_pcm16_wav(&pressure, rate, 5000.0).map_err(|e| e.to_string())?;
    let mut output = std::fs::OpenOptions::new().write(true).create_new(true)
        .open(&args[1]).map_err(|e| e.to_string())?;
    output.write_all(&wav).map_err(|e| e.to_string())?;
    if let Some(path) = args.get(2) {
        let file = std::fs::OpenOptions::new().write(true).create_new(true)
            .open(path).map_err(|e| e.to_string())?;
        let mut csv = std::io::BufWriter::new(file);
        write!(csv, "sample,time_s,mouthpiece_pa,outgoing_pa").map_err(|e| e.to_string())?;
        for mode in 0..modes {
            write!(csv, ",mouth_u{mode}_m3_s").map_err(|e| e.to_string())?;
        }
        writeln!(csv).map_err(|e| e.to_string())?;
        for (sample, ((&p, &a), flow)) in pressure.iter().zip(&outgoing)
            .zip(mouth.chunks(modes)).enumerate()
        {
            write!(csv, "{sample},{:e},{p:e},{a:e}", sample as f64 * dt).map_err(|e| e.to_string())?;
            for q in flow { write!(csv, ",{q:e}").map_err(|e| e.to_string())?; }
            writeln!(csv).map_err(|e| e.to_string())?;
        }
        csv.flush().map_err(|e| e.to_string())?;
    }
    println!("Mouthpiece pressure WAV and optional modal-flow CSV, NOT a microphone; manufactured sources. samples={} rate_hz={} full_scale_pa=5000 clips={} valve_changes={}",
        pressure.len(), rate, clips, sequence.source().voice().lift_log().len());
    for receipt in sequence.receipts() {
        println!("path={} mouth_radius_m={} output_taps={} rungs={}",
            receipt.label, receipt.mouth_radius_m, receipt.n_fft, receipt.rungs.len());
    }
    Ok(())
}
fn main() {
    if let Err(error) = run() { eprintln!("{error}"); std::process::exit(1); }
}
