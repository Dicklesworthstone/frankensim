//! Offline render of the LIVE model from a reference performance CSV
//! (`sample,event,key,value`), for like-for-like comparison with the offline
//! reference binary. Float32 mono WAV, sample = pressure / full scale.
//!
//! usage: render PERFORMANCE.csv SECONDS OUT.wav [--band HZ] [--budget MODES]
//!        [--partials N] [--partial-band HZ] [--full-scale PA] [--scale SUBSET.csv]
//!        [--isolated] [--no-duplex] [--damper-points K] [--trace OUT.csv]
//!        [--reduction MAX,KEEP_LOW,HZ,...] [--max-modes N] [--max-keys N]
//!        [--board FILE.bin]   (a board file written by gen_board, as the browser fetches it)
//! `--scale` + `--isolated` reproduces the reference binary run with the same
//! sparse `--scale` (nothing else attached to the bridge, no retirement).
use fs_live_piano::refmods::geometry;
use fs_live_piano::{Config, Description, LivePiano};
use std::io::Write;

fn main() {
    if let Err(e) = run() {
        eprintln!("render: {e}");
        std::process::exit(1);
    }
}

fn run() -> Result<(), String> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.len() < 3 {
        return Err("usage: render PERFORMANCE.csv SECONDS OUT.wav [options]".into());
    }
    let rate = 48_000u32;
    let mut cfg = Config::live(rate);
    cfg.full_scale_pa = 2.0; // the reference binary's default PCM full scale
    cfg.max_courses = 88;
    cfg.max_dynamic_modes = usize::MAX / 4;
    let mut scale_path: Option<String> = None;
    let mut reduction: Option<String> = None;
    // Default: the embedded board the browser uses. `--band` selects a
    // freshly solved FE slice instead (with `--reduction`, its source band).
    let mut embedded = true;
    let mut trace_path: Option<String> = None;
    let mut board_path: Option<String> = None;
    let mut i = 3;
    while i < args.len() {
        let value = |i: usize| args.get(i + 1).ok_or_else(|| format!("{} needs a value", args[i]));
        let number = |i: usize| value(i).and_then(|v| v.parse::<f64>().map_err(|_| format!("invalid number for {}", args[i])));
        match args[i].as_str() {
            "--band" => { cfg.board_band_hz = number(i)?; embedded = false; i += 1; }
            "--budget" => { cfg.course_mode_budget = number(i)? as usize; i += 1; }
            "--partials" => { cfg.max_partials = number(i)? as usize; i += 1; }
            "--partial-band" => { cfg.partial_band_hz = number(i)?; i += 1; }
            "--full-scale" => { cfg.full_scale_pa = number(i)?; i += 1; }
            "--damper-points" => { cfg.damper_points = number(i)? as usize; i += 1; }
            "--max-modes" => { cfg.max_dynamic_modes = number(i)? as usize; i += 1; }
            "--max-keys" => { cfg.max_courses = number(i)? as usize; i += 1; }
            "--scale" => { scale_path = Some(value(i)?.clone()); i += 1; }
            "--reduction" => { reduction = Some(value(i)?.clone()); i += 1; }
            "--trace" => { trace_path = Some(value(i)?.clone()); i += 1; }
            "--board" => { board_path = Some(value(i)?.clone()); i += 1; }
            "--isolated" => cfg.static_attachment = false,
            "--no-duplex" => cfg.duplex = false,
            other => return Err(format!("unknown option {other}")),
        }
        i += 1;
    }
    let seconds: f64 = args[1].parse().map_err(|_| "invalid duration")?;
    let count = (seconds * f64::from(rate)).round() as usize;
    let text = std::fs::read_to_string(&args[0]).map_err(|e| format!("{}: {e}", args[0]))?;
    let mut events: Vec<(usize, String, u8, f64)> = Vec::new();
    for (line, raw) in text.lines().enumerate() {
        let raw = raw.trim();
        if raw.is_empty() || raw.starts_with('#') || raw.starts_with("sample") { continue; }
        let f: Vec<&str> = raw.split(',').map(str::trim).collect();
        if f.len() != 4 { return Err(format!("line {}: expected sample,event,key,value", line + 1)); }
        events.push((f[0].parse().map_err(|_| "invalid sample")?, f[1].to_owned(),
            f[2].parse().map_err(|_| "invalid key")?, f[3].parse().map_err(|_| "invalid value")?));
    }
    events.sort_by_key(|e| e.0);
    let courses = match &scale_path {
        Some(path) => geometry::read_scale(&std::fs::read_to_string(path).map_err(|e| format!("{path}: {e}"))?)?,
        None => geometry::read_scale(fs_live_piano::kernel::SCALE_CSV)?,
    };
    // `--reduction SPEC` = the offline binary's `--board-reduction SPEC` with
    // `--band` as its `--board-band-hz` source band.
    let desc = match (&reduction, embedded) {
        _ if board_path.is_some() => {
            let path = board_path.as_deref().unwrap_or_default();
            Description::from_board_bytes(courses, &std::fs::read(path).map_err(|e| format!("{path}: {e}"))?)?
        }
        (Some(spec), _) => Description::with_reduction(courses, cfg.board_band_hz, spec)?,
        (None, true) => Description::from_board_bytes(courses, fs_live_piano::kernel::BOARD_BIN)?,
        (None, false) => Description::with_courses(courses, cfg.board_band_hz)?,
    };
    println!("board: {} coordinates; {}", desc.board.len(), desc.board_provenance);
    let mut piano = LivePiano::new(&desc, cfg.clone())?;
    println!("{}", fs_live_piano::describe(&piano));
    let mut out = vec![0.0_f32; count];
    let mut next = 0;
    let scale = 1.0 / cfg.full_scale_pa;
    let trace_key = events.iter().find(|e| e.1 == "note_on").map(|e| e.2);
    let mut trace: Vec<(f64, f64, f64)> = Vec::new();
    let (mut max_modes, mut max_keys) = (0usize, 0usize);
    let start = std::time::Instant::now();
    for (n, sample) in out.iter_mut().enumerate() {
        while next < events.len() && events[next].0 <= n {
            let (_, kind, key, value) = &events[next];
            match kind.as_str() {
                "note_on" => piano.note_on(*key, *value)?,
                "note_off" => piano.note_off(*key),
                "sustain" => piano.set_sustain(*value)?,
                other => return Err(format!("event {other} is not part of the live model")),
            }
            next += 1;
        }
        let pressure = piano.step();
        *sample = (pressure * scale) as f32;
        max_modes = max_modes.max(piano.dynamic_modes());
        max_keys = max_keys.max(piano.sounding_keys());
        if let (Some(key), true) = (trace_key, trace_path.is_some()) {
            trace.push((piano.bridge_velocity(key).unwrap_or(0.0), piano.hammer_force(key).unwrap_or(0.0), pressure));
        }
        if let Some(e) = piano.take_error() { return Err(format!("sample {n}: {e}")); }
    }
    let wall = start.elapsed().as_secs_f64();
    let mut wav = Vec::with_capacity(44 + 4 * count);
    let data = (4 * count) as u32;
    wav.extend_from_slice(b"RIFF");
    wav.extend_from_slice(&(36 + data).to_le_bytes());
    wav.extend_from_slice(b"WAVEfmt ");
    wav.extend_from_slice(&16u32.to_le_bytes());
    wav.extend_from_slice(&3u16.to_le_bytes());
    wav.extend_from_slice(&1u16.to_le_bytes());
    wav.extend_from_slice(&rate.to_le_bytes());
    wav.extend_from_slice(&(rate * 4).to_le_bytes());
    wav.extend_from_slice(&4u16.to_le_bytes());
    wav.extend_from_slice(&32u16.to_le_bytes());
    wav.extend_from_slice(b"data");
    wav.extend_from_slice(&data.to_le_bytes());
    for s in &out { wav.extend_from_slice(&s.to_le_bytes()); }
    std::fs::OpenOptions::new().write(true).create_new(true).open(&args[2])
        .and_then(|mut f| f.write_all(&wav)).map_err(|e| format!("{}: {e}", args[2]))?;
    if let Some(path) = &trace_path {
        let mut text = String::from("sample,bridge_velocity_m_s,hammer_force_n,pressure_pa\n");
        for (n, (v, f, p)) in trace.iter().enumerate() { text.push_str(&format!("{n},{v:.17e},{f:.17e},{p:.17e}\n")); }
        std::fs::OpenOptions::new().write(true).create_new(true).open(path)
            .and_then(|mut f| f.write_all(text.as_bytes())).map_err(|e| format!("{path}: {e}"))?;
    }
    let l = piano.close_ledger();
    let c = piano.counters;
    println!("{seconds} s rendered in {wall:.3} s (real-time factor {:.2}); peak {:.6} Pa; full scale {} Pa; max sounding keys {max_keys}; max dynamic modes {max_modes}.",
        seconds / wall, c.peak_pa, cfg.full_scale_pa);
    println!("Input {:.9} J; stored {:.9} J; dissipated {:.9} J (strings {:.9}, board {:.9}, dampers {:.9}, felt {:.9}, shank {:.9}, catch {:.9}); retired {:.3e} J; residual {:.3e} J.",
        l.input_j, l.stored_j, l.dissipated_j(), l.string_loss_j, l.board_loss_j, l.damper_loss_j, l.felt_loss_j,
        l.shank_loss_j, l.catch_loss_j, l.retired_j, l.residual_j);
    println!("rebuilds {}, retirements {}, forced {}, contact failures {}, clipped samples {}.",
        c.rebuilds, c.retirements, c.forced_retirements, c.contact_failures, c.clipped_samples);
    Ok(())
}
