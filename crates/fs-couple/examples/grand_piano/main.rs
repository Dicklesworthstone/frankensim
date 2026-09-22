//! Geometry-first grand-piano composition of existing FrankenSim owners.
//! `cargo run --release -p fs-couple --example grand_piano -- --preset steinway-d --render piano.wav`
//! Source-derived geometry is approximate, not a measured specimen digital twin.
mod geometry;
mod linear;
mod board;
mod board_geometry;
mod crowned_board;
mod steinway_d;
mod steinway_scale;
mod performance;
use performance::midi;
mod felt;
mod engine;
mod microphone;
mod audio;
mod hammer_materials;

const USAGE: &str = "grand_piano [--render piano.wav] [--scale strings.csv]
    [--preset steinway-d] [--board board.csv | --board-geometry panel.fsb|panel.fss]
    [--hammers materials.fsh] [--dampers estimated | pads.fspd]
    [--concert-pitch 430..450 | --raw-tensions]
    [--mesh-divisions 4..24] [--dump-geometry panel.fsb] [--dump-obj soundboard.obj]
    [--board-band-hz Hz] [--performance events.csv] [--observer-gain Pa/(m^3/s)]
    [--midi performance.mid] [--midi-channel 1..16]
    [--midi-velocity-max-m-s V] [--midi-half-pedal]
    [--microphone x_m,y_m,z_m] [--microphone-right x_m,y_m,z_m] [--diagnostic-volume]
    [--note 21..108] [--velocity m/s] [--duration seconds]
    [--sample-rate Hz] [--substeps 1..16] [--modes 1..512]
    [--dump-scale strings.csv] [--dump-board board.csv]
--preset steinway-d reconstructs the published 17-rib Model D drawing, with
spruce panel, sugar-pine ribs, maple bridges, cut-off bar and 88 bridge stations.
--board-geometry may replace that board while retaining the preset strings,
felt cards and shank mechanics. Its native header selects flat FSB or crowned
3-D shell FSS. An invalid/missing supplied board never falls back to the preset.
Mesh-generation/export controls cannot accompany a supplied board override.
The preset cannot be combined with --board modal CSV.
It uses Chabassier/Durufle's wrapped-string MODEL table (84 notes plus four
estimated extensions), separate per-key hammer force and relaxation cards, and
published shank geometry reduced to rigid rotation plus one bending coordinate.
The shank is not a fitted oscillator: mass, compliance and jack projection come
from its dimensions/material. This is a linearized reduction, not the full action.
Lengths, effective winding mass and EI remain fixed. Preset tensions are tuned
to A4=440 Hz by default to compensate rounded source values; --raw-tensions
preserves the published table. --concert-pitch tunes first partials by changing
physical tension, not oscillator frequencies. This is not a stretch-tuning fit.
An explicit --scale always supplies the geometry/masses; absent --concert-pitch,
its tensions are preserved even with a preset. Preset hammer voicing still applies
unless --hammers supplies a complete per-key WoolFelt/Prony material file.
--hammers requires --render and a card for every key in the supplied scale, not
only the struck keys. Missing, duplicate or invalid cards refuse without fallback.
These cards change physical contact forces and relaxation, not an output EQ.
The preset shank and the scale's hammer mass/patch geometry remain unchanged.
See HAMMERS.md for the SI format; importing values does not certify measurements.
--dampers selects finite-footprint viscous pads instead of the default point
damper. 'estimated' declares approximate spans and drag; a file must cover
every scale key with a pad or explicit free row. It requires --render and
uses existing MIDI/CSV key, sustain and sostenuto controls. See DAMPERS.md.
This is spatial drag, not falling-pad or hysteretic felt contact mechanics.
--dump-geometry/--dump-obj export the board model; export alone skips eigenanalysis.
Thickness taper, material constants and key assignment include explicit estimates.
Default per-key WoolFelt loading envelopes are source-derived; crush/unloading
parameters, felt patch geometry and Prony time constants remain estimates.
--note performs a single-key study; otherwise the demo also plays a chord of
available keys. --velocity overrides the three demo hammer launch speeds.
Velocity is POST-ESCAPEMENT hammer velocity, not MIDI velocity or key motion.
Geometric boards assemble a flat orthotropic plate or a supplied crowned shell.
Crowned structures retain 3-D motion, grain and eccentric ribs/bridges; their
radiation is a projected flat-baffle approximation, not 3-D exterior BEM.
Crown is reference geometry, not solved downbearing; see CROWNED_BOARD.md.
The explicit frequency band admits at most 128 modes; --modes admits up to 512
string partials. Neither budget nor mesh refinement alone is a convergence or
real-time claim.
--performance uses sample,event,key,value CSV instead of the demo and cannot be
combined with --note or --velocity. note_on values are hammer velocity in m/s.
With the preset, jack_staccato and jack_legato instead take peak force in N at
the physical jack station, with 7/100 ms pulses and 1.5 mm let-off. For example:
0,jack_staccato,27,70
48000,jack_legato,69,30
Use note_off events to release keys before repeating them. Shank damping and
the backcheck remain estimates; jack timing is resolved at the mechanical rate.
--midi replaces the demo/--performance with a format-0/1 Standard MIDI File.
It selects channel 1 by default. Velocity 127 maps to 4.5 m/s by default;
--midi-velocity-max-m-s explicitly changes this linear, uncalibrated hammer
launch mapping, not audio gain. CC64/66/67 drive the existing three pedals;
--midi-half-pedal opts into CC64/127 travel instead of the default switch.
The complete score must fit --duration, leaving room for acoustic ringdown.
Unsupported physical controls refuse or are reported as ignored; see MIDI.md.
Geometric boards default to a spatial Rayleigh half-space pressure microphone
at (0.675,1,1) metres in the mesh coordinate system. --microphone moves it.
--microphone-right adds a second spatial receiver and writes left/right stereo
from one mechanical performance. The original/default microphone is left;
no panning, duplicated mono, peak normalization or second physics run. See STEREO.md.
This assumes an infinite baffle, with no lid/room scattering or air backreaction.
Pressure uses every mechanics substep and causal anti-alias filtering before
output-rate propagation. At 4x oversampling the filter adds 44 audio samples of
latency, in addition to acoustic travel time. Histories persist across blocks.
--diagnostic-volume retains the old volume-velocity observer; --observer-gain
applies only to that diagnostic, not to physical microphone pressure.";

#[derive(Debug)]
struct Options {
    render: Option<String>, scale: Option<String>, board: Option<String>,
    board_geometry: Option<String>, performance: Option<String>, preset: Option<String>,
    hammers: Option<String>, dampers: Option<String>,
    midi: Option<String>, midi_mapping: midi::Mapping,
    concert_pitch: Option<f64>, raw_tensions: bool,
    mesh_divisions: usize, dump_geometry: Option<String>, dump_obj: Option<String>,
    board_band_hz: f64, observer_gain: f64,
    microphone: Option<[f64; 3]>, microphone_right: Option<[f64; 3]>, diagnostic_volume: bool,
    dump_scale: Option<String>, dump_board: Option<String>,
    note: Option<u8>, velocity: Option<f64>, duration: f64,
    sample_rate: u32, substeps: usize, modes: usize, help: bool,
}
impl Default for Options {
    fn default() -> Self {
        Self { render: None, scale: None, board: None, board_geometry: None,
            performance: None, preset: None, hammers: None, dampers: None, concert_pitch: None, raw_tensions: false,
            midi: None, midi_mapping: midi::Mapping::default(),
            mesh_divisions: 8, dump_geometry: None, dump_obj: None,
            board_band_hz: 400.0, observer_gain: 10_000.0, dump_scale: None,
            microphone: None, microphone_right: None, diagnostic_volume: false,
            dump_board: None, note: None, velocity: None, duration: 6.0,
            sample_rate: 48_000, substeps: 4, modes: 24, help: false }
    }
}
impl Options {
    fn parse(args: &[String]) -> Result<Self, String> {
        let mut options = Self::default();
        let mut seen = std::collections::BTreeSet::new();
        let mut args = args.iter();
        while let Some(flag) = args.next() {
            if flag == "--help" || flag == "-h" { options.help = true; continue; }
            if !seen.insert(flag.as_str()) { return Err(format!("duplicate option {flag}")); }
            if flag == "--diagnostic-volume" { options.diagnostic_volume = true; continue; }
            if flag == "--raw-tensions" { options.raw_tensions = true; continue; }
            if flag == "--midi-half-pedal" { options.midi_mapping.continuous_sustain = true; continue; }
            let value = args.next().ok_or_else(|| format!("missing value for {flag}"))?;
            let invalid = || format!("invalid value for {flag}: {value}");
            match flag.as_str() {
                "--render" => options.render = Some(value.clone()),
                "--scale" => options.scale = Some(value.clone()),
                "--board" => options.board = Some(value.clone()),
                "--board-geometry" => options.board_geometry = Some(value.clone()),
                "--preset" => options.preset = Some(value.clone()),
                "--hammers" => options.hammers = Some(value.clone()),
                "--dampers" => options.dampers = Some(value.clone()),
                "--concert-pitch" => options.concert_pitch = Some(value.parse().map_err(|_| invalid())?),
                "--mesh-divisions" => options.mesh_divisions = value.parse().map_err(|_| invalid())?,
                "--dump-geometry" => options.dump_geometry = Some(value.clone()),
                "--dump-obj" => options.dump_obj = Some(value.clone()),
                "--performance" => options.performance = Some(value.clone()),
                "--midi" => options.midi = Some(value.clone()),
                "--midi-channel" => {
                    options.midi_mapping.channel = value.parse::<u8>().map_err(|_| invalid())?
                        .checked_sub(1).ok_or_else(invalid)?;
                }
                "--midi-velocity-max-m-s" => options.midi_mapping.maximum_velocity_m_s =
                    value.parse().map_err(|_| invalid())?,
                "--microphone" | "--microphone-right" => {
                    let values = value.split(',').map(str::parse::<f64>).collect::<Result<Vec<_>,_>>()
                        .map_err(|_| invalid())?;
                    if values.len() != 3 { return Err(invalid()); }
                    let position=Some([values[0], values[1], values[2]]);
                    if flag=="--microphone" {options.microphone=position;} else {options.microphone_right=position;}
                }
                "--board-band-hz" => options.board_band_hz = value.parse().map_err(|_| invalid())?,
                "--observer-gain" => options.observer_gain = value.parse().map_err(|_| invalid())?,
                "--dump-scale" => options.dump_scale = Some(value.clone()),
                "--dump-board" => options.dump_board = Some(value.clone()),
                "--note" => options.note = Some(value.parse().map_err(|_| invalid())?),
                "--velocity" => options.velocity = Some(value.parse().map_err(|_| invalid())?),
                "--duration" => options.duration = value.parse().map_err(|_| invalid())?,
                "--sample-rate" => options.sample_rate = value.parse().map_err(|_| invalid())?,
                "--substeps" => options.substeps = value.parse().map_err(|_| invalid())?,
                "--modes" => options.modes = value.parse().map_err(|_| invalid())?,
                _ => return Err(format!("unknown option {flag}\n{USAGE}")),
            }
        }
        if options.note.is_some_and(|n| !(21..=108).contains(&n))
            || options.velocity.is_some_and(|v| !v.is_finite() || v <= 0.0 || v > 8.0)
            || !options.duration.is_finite() || !(0.001..=120.0).contains(&options.duration)
            || !(8_000..=192_000).contains(&options.sample_rate)
            || !(1..=16).contains(&options.substeps) || !(1..=linear::MAX_STRING_MODES).contains(&options.modes) {
            return Err("render control outside its finite admitted range".into());
        }
        if options.concert_pitch.is_some_and(|f| !f.is_finite() || !(430.0..=450.0).contains(&f))
            || (options.raw_tensions && (options.preset.is_none() || options.concert_pitch.is_some())) {
            return Err("concert pitch must be 430..450 Hz; --raw-tensions requires a preset and excludes --concert-pitch".into());
        }
        if options.board.is_some() && (options.board_geometry.is_some() || options.preset.is_some()) {
            return Err("modal board CSV excludes a preset and geometric board; --board-geometry may override a preset board".into());
        }
        if options.preset.as_deref().is_some_and(|p| p != "steinway-d") {
            return Err("unknown piano preset; available: steinway-d".into());
        }
        if options.hammers.is_some() && options.render.is_none() {
            return Err("--hammers requires --render; material input is not an export-only option".into());
        }
        if options.dampers.as_ref().is_some_and(|s| s.trim().is_empty() || options.render.is_none()) {
            return Err("--dampers requires --render and either estimated or a nonempty specification path".into());
        }
        if !(4..=24).contains(&options.mesh_divisions)
            || (!options.uses_preset_board() && (seen.contains("--mesh-divisions")
                || options.dump_geometry.is_some() || options.dump_obj.is_some())) {
            return Err("mesh/export controls require --preset steinway-d without a supplied board override; divisions must be 4..24".into());
        }
        if seen.contains("--board-band-hz") && options.board_geometry.is_none() && options.preset.is_none() {
            return Err("--board-band-hz requires --board-geometry or --preset".into());
        }
        if !options.board_band_hz.is_finite() || options.board_band_hz <= 0.0
            || options.board_band_hz >= 0.45 * f64::from(options.sample_rate)
            || !options.observer_gain.is_finite() || options.observer_gain <= 0.0 {
            return Err("invalid board frequency band or diagnostic observer gain".into());
        }
        if options.performance.is_some()
            && (options.render.is_none() || options.note.is_some() || options.velocity.is_some()) {
            return Err("--performance requires --render and replaces --note/--velocity demo controls".into());
        }
        if options.midi.is_some() && (options.render.is_none() || options.performance.is_some()
            || options.note.is_some() || options.velocity.is_some()) {
            return Err("--midi requires --render and excludes --performance/--note/--velocity".into());
        }
        let mapping = options.midi_mapping;
        if mapping.channel > 15 || !mapping.maximum_velocity_m_s.is_finite()
            || mapping.maximum_velocity_m_s <= 0.0 || mapping.maximum_velocity_m_s > 8.0
            || mapping.maximum_velocity_m_s / 127.0 == 0.0
            || (options.midi.is_none() && ["--midi-channel", "--midi-velocity-max-m-s", "--midi-half-pedal"]
                .iter().any(|flag| seen.contains(flag))) {
            return Err("MIDI controls require --midi, channel 1..16 and finite maximum hammer velocity in (0,8] m/s".into());
        }
        let geometric = options.preset.is_some() || options.board_geometry.is_some();
        if [options.microphone,options.microphone_right].iter().flatten()
            .any(|p| p.iter().any(|x|!x.is_finite()) || p[2]<0.05)
            || ((options.microphone.is_some() || options.microphone_right.is_some())
                && (!geometric || options.render.is_none() || options.diagnostic_volume)) {
            return Err("microphones need a geometric render, finite x,y,z with z>=0.05, and no --diagnostic-volume".into());
        }
        if geometric && !options.diagnostic_volume && seen.contains("--observer-gain") {
            return Err("--observer-gain requires --diagnostic-volume for a geometric board".into());
        }
        // Do not overwrite the very measurements that a render was asked to use.
        let inputs = [options.scale.as_ref(), options.board.as_ref(),
            options.board_geometry.as_ref(), options.performance.as_ref(), options.hammers.as_ref(), options.midi.as_ref(),
            options.dampers.as_ref().filter(|s| s.as_str() != "estimated")];
        let outputs = [options.render.as_ref(), options.dump_scale.as_ref(), options.dump_board.as_ref(),
            options.dump_geometry.as_ref(), options.dump_obj.as_ref()];
        for (i, output) in outputs.iter().enumerate() {
            if let Some(path) = output {
                if path.is_empty() || inputs.iter().flatten().any(|input| input == path)
                    || outputs[..i].iter().flatten().any(|previous| previous == path) {
                    return Err("output paths must be distinct from inputs and each other".into());
                }
            }
        }
        Ok(options)
    }
    fn tuning_hz(&self) -> Option<f64> {
        self.concert_pitch.or_else(||
            (self.preset.is_some() && self.scale.is_none() && !self.raw_tensions).then_some(440.0))
    }
    fn uses_preset_board(&self) -> bool {
        self.preset.is_some() && self.board_geometry.is_none()
    }
}

fn load_scale(text: Option<&str>) -> Result<Vec<geometry::Course>, String> {
    match text { Some(text) => geometry::read_scale(text), None => geometry::demonstration_scale() }
}
fn selected_scale(text: Option<&str>, options: &Options) -> Result<Vec<geometry::Course>, String> {
    let mut scale = if text.is_none() && options.preset.is_some() {
        steinway_scale::courses()?
    } else { load_scale(text)? };
    if let Some(reference) = options.tuning_hz() {
        for c in &mut scale {
            let target = reference * fs_math::det::pow(2.0, (f64::from(c.midi) - 69.0) / 12.0);
            let cents = 1200.0 * fs_math::det::ln(target / c.partial_hz(1, c.tension_n))
                / std::f64::consts::LN_2;
            c.tension_n = c.tension_at_cents(cents)
                .map_err(|e| format!("key {} tension retuning: {e}", c.midi))?;
        }
    }
    Ok(scale)
}
fn prepare_instrument(scale: Vec<geometry::Course>, modes: &[linear::BoardMode],
    options: &Options) -> Result<engine::Instrument, String> {
    let dampers = match options.dampers.as_deref() {
        None => None,
        Some("estimated") => Some(linear::dampers::Specification::estimated(&scale)?),
        Some(path) => Some(linear::dampers::Specification::load(path, &scale)?),
    };
    let text = options.hammers.as_ref().map(|path| std::fs::read_to_string(path)
        .map_err(|e| format!("{path}: {e}"))).transpose()?;
    let mut piano = prepare_instrument_with_hammers(scale, modes, options, text.as_deref())?;
    if let Some(spec) = &dampers { piano.configure_dampers(spec)?; }
    Ok(piano)
}
fn prepare_instrument_with_hammers(scale: Vec<geometry::Course>, modes: &[linear::BoardMode],
    options: &Options, text: Option<&str>) -> Result<engine::Instrument, String> {
    let imported = text.map(|text| hammer_materials::read(text,
        &scale.iter().map(|c| c.midi).collect::<Vec<_>>())).transpose()?;
    if options.preset.is_some() {
        let materials = match imported {
            Some(materials) => materials,
            None => scale.iter().map(steinway_scale::hammer_material).collect::<Result<Vec<_>,_>>()?,
        };
        engine::Instrument::new_with_course_shanks(scale, modes, options.sample_rate,
            options.substeps, options.modes, true, materials, engine::ShankGeometry::published())
    } else if let Some(materials) = imported {
        engine::Instrument::new_with_course_felts(scale, modes, options.sample_rate,
            options.substeps, options.modes, true, materials)
    } else {
        engine::Instrument::new(scale, modes, options.sample_rate, options.substeps, options.modes, true)
    }
}
fn load_board(text: Option<&str>, scale: &[geometry::Course]) -> Result<Vec<linear::BoardMode>, String> {
    match text {
        Some(text) => board::read(text, &scale.iter().map(|c| c.midi).collect::<Vec<_>>()),
        None => Ok(board::demonstration()),
    }
}
fn prepare_geometric_board(text: &str, keys: &[u8], band_hz: f64)
    -> Result<board_geometry::PreparedBoard, String> {
    if crowned_board::is_crowned(text) {
        crowned_board::CrownedBoard::read(text)?.prepare(keys, band_hz)
    } else { board_geometry::BoardGeometry::read(text)?.prepare(keys, band_hz) }
}
/// Export only the admitted keys. An absent measurement must not become an
/// apparently measured zero bridge coefficient when the table is re-imported.
fn write_board_for_scale(modes: &[linear::BoardMode], scale: &[geometry::Course]) -> String {
    let all = board::write(modes);
    let mut out = String::new();
    for row in all.lines() {
        let admitted = row == board::HEADER || row.split(',').nth(4)
            .and_then(|key| key.parse::<u8>().ok())
            .is_some_and(|key| scale.iter().any(|c| c.midi == key));
        if admitted { out.push_str(row); out.push('\n'); }
    }
    out
}
fn study_key(scale: &[geometry::Course], requested: Option<u8>) -> Result<u8, String> {
    if let Some(key) = requested {
        return scale.iter().any(|c| c.midi == key).then_some(key)
            .ok_or_else(|| format!("requested key {key} is absent from the input scale"));
    }
    scale.iter().min_by_key(|c| c.midi.abs_diff(69)).map(|c| c.midi)
        .ok_or_else(|| "empty input scale".into())
}

fn render(path: &str, scale: Vec<geometry::Course>, modes: &[linear::BoardMode],
    surface: Option<&[board_geometry::SurfaceSample]>, options: &Options) -> Result<(), String> {
    study_key(&scale, options.note)?;
    let keys: Vec<u8> = scale.iter().map(|c| c.midi).collect();
    let rate = options.sample_rate;
    let count = (options.duration * f64::from(rate)).round() as u32;
    let score = if let Some(path) = &options.midi {
        let parsed = midi::load(path, &keys, rate, u64::from(count), options.midi_mapping)?;
        let report = &parsed.report;
        println!("MIDI: {} tracks, channel {}, {} hammer launches, end sample {}, {} end releases; {} other-channel messages, {} unsupported channel messages and {} SysEx events ignored.",
            report.tracks, options.midi_mapping.channel + 1, report.selected_note_ons,
            report.end_sample, report.end_releases, report.other_channel_messages,
            report.ignored_channel_messages, report.ignored_sysex_events);
        println!("Uncalibrated MIDI mapping: velocity 127 -> {} m/s, linear hammer launch; sustain mapping {}. No output gain, sample bank or pitch-wheel substitution.",
            options.midi_mapping.maximum_velocity_m_s,
            if options.midi_mapping.continuous_sustain { "CC64/127 travel" } else { "switch at 64" });
        performance::Performance::from_events(parsed.events, &keys, u64::from(count))?
    } else { match &options.performance {
        Some(path) => performance::Performance::read(
            &std::fs::read_to_string(path).map_err(|e| format!("{path}: {e}"))?,
            &keys, u64::from(count))?,
        None => performance::Performance::demonstration(&keys, rate, u64::from(count),
            options.note, options.velocity)?,
    }};
    let piano = prepare_instrument(scale, modes, options)?;
    debug_assert_eq!(piano.sample_rate(), rate);
    if let Some((strings,cells)) = piano.damper_resolution() {
        println!("Spatial viscous dampers: {} speaking-string pads, {} quadrature stations; source {}. No measured pad/action or real-time claim; see DAMPERS.md.",
            strings, cells, options.dampers.as_deref().unwrap_or("supplied"));
    }
    let surface = if options.diagnostic_volume { None } else { surface };
    let left=options.microphone.unwrap_or([0.675,1.0,1.0]);
    let mut stream = match options.microphone_right {
        Some(right)=>audio::AudioStream::new_stereo(piano,score,
            surface.ok_or("stereo pressure requires a geometric soundboard surface")?,
            [left,right],fs_bem::helmholtz::Medium::air())?,
        None=>audio::AudioStream::new(piano,score,surface,left,
            fs_bem::helmholtz::Medium::air(),options.observer_gain)?,
    };
    if let Some(mic) = stream.microphone() {
        println!("Rayleigh pressure microphone at {:?} m; propagation {:?} samples plus {} anti-alias delay samples; {} modal multiplies/output sample.",
            mic.position_m, mic.delay_samples, mic.filter_delay_samples(), mic.multiply_adds_per_sample());
    }
    if let Some([_,right])=stream.stereo_microphones() {
        println!("Right Rayleigh microphone at {:?} m; propagation {:?} samples plus {} anti-alias delay samples. Independent spatial pressure; shared mechanics.",
            right.position_m,right.delay_samples,right.filter_delay_samples());
    }
    let channels=stream.channels();
    let mut pressure = vec![0.0; count as usize*channels];
    let start = std::time::Instant::now();
    for block in pressure.chunks_mut(256*channels) {
        stream.render_interleaved_block(block).map_err(|e| e.to_string())?;
    }
    let elapsed = start.elapsed().as_secs_f64();
    let peak = pressure.iter().fold(0.0_f64, |a, p| a.max(p.abs()));
    let seconds = f64::from(count) / f64::from(rate);
    let (wav, clips) = fs_couple::pcm_wav::encode_pcm16_wav_interleaved(&pressure, rate, channels as u16, 2.0).map_err(|e| e.to_string())?;
    std::fs::write(path, wav).map_err(|e| format!("{path}: {e}"))?;
    if stream.microphone().is_some() {
        println!("Computed half-space pressure in Pa; PCM full scale 2 Pa, no peak normalization. Infinite baffle; no room/lid scattering, radiation loading or measured-SPL calibration.");
    } else {
        println!("Diagnostic volume-velocity observer, gain {} Pa/(m^3/s); no peak normalization or calibrated SPL claim.", options.observer_gain);
    }
    let piano = stream.instrument();
    println!("{} string modes; {} board modes; {} above-band duplex segments omitted from dynamic retention (static attachment retained).",
        piano.bank.modes.len(), piano.bank.board_count, piano.bank.omitted_duplex_modes);
    println!("{seconds:.6} s, {channels} channel(s) rendered in {elapsed:.6} s; wall/audio ratio {:.4}; peak {peak:.6} Pa-equivalent; {clips} PCM clips.", elapsed / seconds);
    println!("Input {:.9} J; stored {:.9} J; component losses {:.9} J; closure {:.3e} J; worst substep defect {:.3e} J.",
        piano.accounting.input_work_j, piano.energy_j(), piano.accounting.dissipated_j(),
        piano.accounting.input_work_j - piano.energy_j() - piano.accounting.dissipated_j(), piano.accounting.max_balance_error_j);
    println!("Felt loss {:.9} J, including {:.9} J time-dependent relaxation; shank damping {:.9} J.",
        piano.accounting.felt_loss_j, piano.accounting.felt_relaxation_loss_j, piano.accounting.shank_loss_j);
    if piano.damper_resolution().is_some() {
        println!("Spatial damper loss {:.9} J, already included in component losses; no output-envelope damping.",
            piano.accounting.damper_loss_j);
    }
    Ok(())
}

fn run() -> Result<(), String> {
    let options = Options::parse(&std::env::args().skip(1).collect::<Vec<_>>())?;
    if options.help { println!("{USAGE}"); return Ok(()); }
    let read = |path: &String| std::fs::read_to_string(path).map_err(|e| format!("{path}: {e}"));
    let scale_text = options.scale.as_ref().map(read).transpose()?;
    let board_text = options.board.as_ref().map(read).transpose()?;
    let scale = selected_scale(scale_text.as_deref(), &options)?;
    let scale_source = options.scale.as_deref().unwrap_or(if options.preset.is_some() {
        "RT-0425 Appendix A wrapped-string MODEL: 84 published courses plus four estimated extensions"
    } else { "ESTIMATED demonstration" });
    let tuning_source = options.tuning_hz().map_or_else(|| "input tensions preserved".to_owned(),
        |f| format!("tensions adjusted to A4={f} Hz first-partial equal temperament; L, mass and EI preserved"));
    let preset = options.uses_preset_board().then(|| steinway_d::build(options.mesh_divisions)).transpose()?;
    if let Some(preset) = &preset {
        if let Some(path) = &options.dump_geometry { std::fs::write(path, &preset.geometry).map_err(|e| format!("{path}: {e}"))?; }
        if let Some(path) = &options.dump_obj { std::fs::write(path, &preset.obj).map_err(|e| format!("{path}: {e}"))?; }
        if options.render.is_none() && options.dump_board.is_none() && options.dump_scale.is_none()
            && (options.dump_geometry.is_some() || options.dump_obj.is_some()) {
            println!("Exported source-derived Model D: tapered panel, 17 ribs, maple bridges, cut-off bar and 88 bridge stations. No eigenanalysis was needed.");
            return Ok(());
        }
    }
    let geometry_text = match (&preset, &options.board_geometry) {
        (_, Some(path)) => Some(read(path)?),
        (Some(p), _) => Some(p.geometry.clone()),
        _ => None,
    };
    let (modes, board_source, surface) = if let Some(text) = &geometry_text {
        let start = std::time::Instant::now();
        let prepared = prepare_geometric_board(text, &scale.iter().map(|c| c.midi).collect::<Vec<_>>(),
            options.board_band_hz)?;
        let model_name = if crowned_board::is_crowned(text) { "Crowned shell" } else { "Flat plate" };
        println!("{model_name} {:.6} m^2, {:.6} kg (panel+ribs/bridges), {} free DOFs, {} modes in (0,{}] Hz; preparation {:.6} s.",
            prepared.area_m2, prepared.mass_kg, prepared.free_dofs, prepared.modes.len(),
            options.board_band_hz, start.elapsed().as_secs_f64());
        for (i, interval) in prepared.frequency_intervals_hz.iter().enumerate() {
            println!("board mode {i}: [{:.9}, {:.9}] Hz", interval.0, interval.1);
        }
        (prepared.modes, format!("GEOMETRY-DERIVED {model_name}; {}; rim compliance not modeled", prepared.provenance),
            Some(prepared.surface))
    } else {
        (load_board(board_text.as_deref(), &scale)?, options.board.as_deref()
            .unwrap_or("AUTHORED illustrative modes; not measured Steinway geometry").to_owned(), None)
    };
    if modes.iter().any(|m| m.frequency_hz >= 0.45 * f64::from(options.sample_rate)) {
        return Err("soundboard mode at/above output retention ceiling; use an explicitly reduced board".into());
    }
    study_key(&scale, options.note)?;
    println!("String scale: {scale_source}; {tuning_source}.");
    println!("Soundboard: {board_source}.");
    if let Some(path) = &options.hammers {
        println!("Hammer materials supplied by {path}; no automatic coupon-fit or calibration claim. Scale mass/patch geometry unchanged.");
    } else if options.preset.is_some() {
        println!("Per-key source-derived hammer loading envelopes; estimated unloading/crush and tangent-scaled Prony relaxation.");
    }
    if options.preset.is_some() {
        println!("Published shank geometry -> rigid rotation + bending; reciprocal jack port and 1.5 mm let-off. Linearized action fragment; damping/backcheck estimated.");
    }
    println!("Source authority belongs to the inputs, not the model name; imported files are not independently certified measurements.");
    if let Some(path) = &options.dump_scale {
        std::fs::write(path, format!("# Source: {scale_source}; {tuning_source}\n{}", geometry::write_scale(&scale))).map_err(|e| e.to_string())?;
    }
    if let Some(path) = &options.dump_board {
        std::fs::write(path, format!("# Source: {board_source}\n{}", write_board_for_scale(&modes, &scale)))
            .map_err(|e| e.to_string())?;
    }
    if let Some(path) = &options.render { return render(path, scale, &modes, surface.as_deref(), &options); }
    println!("Model D published envelope: {} x {} m; board {} -> {} m (center -> edge).",
        geometry::D_LENGTH_M, geometry::D_WIDTH_M, geometry::D_BOARD_CENTER_M, geometry::D_BOARD_EDGE_M);
    println!("{} courses; {} speaking strings; {} board modes.", scale.len(), scale.iter().map(|c| c.unison).sum::<usize>(), modes.len());
    for c in scale { println!("key {}: L={:.4} m, mu={:.6} kg/m, T={:.2} N, f1={:.3} Hz", c.midi, c.length_m, c.linear_density_kg_m, c.tension_n, c.partial_hz(1, c.tension_n)); }
    Ok(())
}
fn main() {
    if let Err(error) = run() { eprintln!("grand_piano: {error}"); std::process::exit(1); }
}

#[cfg(test)]
#[path = "crowned_render_tests.rs"]
mod crowned_render_tests;

#[cfg(test)]
mod render_tests {
    use super::*;
    fn options(args: &[&str]) -> Result<Options, String> {
        Options::parse(&args.iter().map(|s| (*s).to_owned()).collect::<Vec<_>>())
    }
    #[test]
    fn damper_selection_composes_with_physical_and_midi_inputs_and_protects_its_source() {
        let o = options(&["--preset", "steinway-d", "--midi", "score.mid", "--render", "piano.wav",
            "--dampers", "pads.fspd", "--hammers", "felt.fsh"]).unwrap();
        assert_eq!(o.dampers.as_deref(), Some("pads.fspd"));
        assert!(options(&["--render", "piano.wav", "--dampers", "estimated"]).is_ok());
        assert_eq!(options(&[]).unwrap().dampers, None);
        for args in [vec!["--dampers"], vec!["--dampers", "estimated"],
            vec!["--render", "piano.wav", "--dampers", ""],
            vec!["--render", "pads.fspd", "--dampers", "pads.fspd"],
            vec!["--render", "piano.wav", "--dampers", "pads.fspd", "--dump-scale", "pads.fspd"],
            vec!["--render", "piano.wav", "--dampers", "estimated", "--dampers", "other.fspd"]] {
            assert!(options(&args).is_err(), "accepted {args:?}");
        }
    }
    #[test]
    fn stereo_receiver_is_explicit_physical_input_and_preserves_existing_controls() {
        let o=options(&["--preset","steinway-d","--render","stereo.wav",
            "--microphone","0.2,0.8,1","--microphone-right","1.2,0.8,1",
            "--midi","score.mid","--dampers","estimated"]).unwrap();
        assert_eq!(o.microphone,Some([0.2,0.8,1.0]));
        assert_eq!(o.microphone_right,Some([1.2,0.8,1.0]));
        assert_eq!(o.midi.as_deref(),Some("score.mid"));
        assert!(options(&["--preset","steinway-d","--render","s.wav",
            "--microphone-right","1,1,1"]).unwrap().microphone.is_none());
        for args in [vec!["--microphone-right","1,1,1"],
            vec!["--render","s.wav","--microphone-right","1,1,1"],
            vec!["--preset","steinway-d","--render","s.wav","--microphone-right","1,1,NaN"],
            vec!["--preset","steinway-d","--render","s.wav","--microphone-right","1,1,0"],
            vec!["--preset","steinway-d","--render","s.wav","--microphone-right","1,1"],
            vec!["--preset","steinway-d","--render","s.wav","--microphone-right","1,1,1","--diagnostic-volume"],
            vec!["--preset","steinway-d","--render","s.wav","--microphone-right","1,1,1","--microphone-right","2,1,1"]] {
            assert!(options(&args).is_err(),"accepted {args:?}");
        }
    }
    #[test]
    fn midi_options_compose_with_physical_inputs_without_overwriting_the_score() {
        let o = options(&["--preset", "steinway-d", "--midi", "score.mid", "--render", "piano.wav",
            "--midi-channel", "16", "--midi-velocity-max-m-s", "2", "--midi-half-pedal",
            "--hammers", "felt.fsh"]).unwrap();
        assert_eq!(o.midi.as_deref(), Some("score.mid"));
        assert_eq!(o.midi_mapping.channel, 15);
        assert_eq!(o.midi_mapping.maximum_velocity_m_s, 2.0);
        assert!(o.midi_mapping.continuous_sustain);
        for args in [vec!["--midi", "score.mid"], vec!["--midi-channel", "1"],
            vec!["--midi-half-pedal"], vec!["--midi-velocity-max-m-s", "2"],
            vec!["--midi", "score.mid", "--render", "score.mid"],
            vec!["--midi", "score.mid", "--render", "p.wav", "--performance", "p.csv"],
            vec!["--midi", "score.mid", "--render", "p.wav", "--note", "69"],
            vec!["--midi", "score.mid", "--render", "p.wav", "--velocity", "2"],
            vec!["--midi", "score.mid", "--render", "p.wav", "--midi-channel", "0"],
            vec!["--midi", "score.mid", "--render", "p.wav", "--midi-channel", "17"],
            vec!["--midi", "score.mid", "--render", "p.wav", "--midi-velocity-max-m-s", "NaN"],
            vec!["--midi", "score.mid", "--render", "p.wav", "--midi-velocity-max-m-s", "9"]] {
            assert!(options(&args).is_err(), "accepted {args:?}");
        }
    }
    #[test]
    fn midi_and_si_csv_drive_identical_preset_contact_pedals_and_energy() {
        // PPQN 480, default tempo: four ticks -> sample 200 at 48 kHz.
        // End-of-track at tick 12 releases sustain at sample 600, not PCM.
        let bytes = b"MThd\0\0\0\x06\0\0\0\x01\x01\xe0MTrk\0\0\0\x10\
            \0\xb0\x40\x7f\0\x90\x45\x7f\x04\x90\x45\0\x08\xff\x2f\0";
        let mut o = options(&["--preset", "steinway-d", "--render", "piano.wav"]).unwrap(); o.modes = 12;
        for spatial in [false,true] {
        o.dampers = spatial.then(|| "estimated".to_owned());
        let course = selected_scale(None, &o).unwrap()[48];
        let mut a = prepare_instrument(vec![course], &board::demonstration(), &o).unwrap();
        let mut b = prepare_instrument(vec![course], &board::demonstration(), &o).unwrap();
        assert_eq!(a.damper_resolution().is_some(), spatial);
        let decoded = midi::read(bytes, &[69], 48_000, 1500,
            midi::Mapping { maximum_velocity_m_s: 2.0, ..midi::Mapping::default() }).unwrap();
        assert_eq!(decoded.report.end_releases, 1);
        let mut imported = performance::Performance::from_events(decoded.events, &[69], 1500).unwrap();
        let mut csv = performance::Performance::read("sample,event,key,value\n0,sustain,0,1\n\
            0,note_on,69,2\n200,note_off,69,0\n600,sustain,0,0\n", &[69], 1500).unwrap();
        for sample in 0..1500 {
            imported.dispatch(sample, &mut a).unwrap(); csv.dispatch(sample, &mut b).unwrap();
            assert_eq!(a.step().unwrap().to_bits(), b.step().unwrap().to_bits());
        }
        assert!(a.accounting.felt_loss_j > 0.0);
        assert_eq!(a.accounting.input_work_j.to_bits(), b.accounting.input_work_j.to_bits());
        assert!((a.accounting.input_work_j - a.energy_j() - a.accounting.dissipated_j()).abs() < 1e-7);
        assert!(a.accounting.damper_loss_j > 0.0);
        }
    }
    #[test]
    fn rendering_and_both_physical_imports_are_composable() {
        let o = options(&["--scale", "strings.csv", "--render", "piano.wav", "--board", "board.csv",
            "--sample-rate", "44100", "--substeps", "2", "--modes", "40", "--note", "21"]).unwrap();
        assert_eq!(o.render.as_deref(), Some("piano.wav"));
        assert_eq!(o.scale.as_deref(), Some("strings.csv"));
        assert_eq!(o.board.as_deref(), Some("board.csv"));
        assert_eq!((o.sample_rate, o.substeps, o.modes, o.note), (44_100, 2, 40, Some(21)));
    }
    #[test]
    fn imported_subset_reaches_the_prepared_instrument_without_substitution() {
        let mut course = geometry::demonstration_scale().unwrap()[0];
        course.length_m *= 0.97;
        let scale = load_scale(Some(&geometry::write_scale(&[course]))).unwrap();
        assert_eq!(scale, vec![course]);
        assert_eq!(study_key(&scale, None).unwrap(), 21);
        assert!(study_key(&scale, Some(69)).is_err());
        let mut input = board::demonstration();
        input.truncate(1);
        input[0].frequency_hz = 123.0;
        let modes = load_board(Some(&board::write(&input)), &scale).unwrap();
        assert_eq!(modes[0].frequency_hz, 123.0);
        let mut piano = engine::Instrument::new(scale, &modes, 8_000, 4, 4, true).unwrap();
        assert_eq!(piano.bank.board_count, 1);
        piano.note_on(21, 1.0).unwrap();
        for _ in 0..16 { assert!(piano.step().unwrap().is_finite()); }
    }
    #[test]
    fn invalid_inputs_never_fall_back_to_estimates() {
        assert!(load_scale(Some("bad CSV")).is_err());
        let scale = geometry::demonstration_scale().unwrap();
        assert!(load_board(Some(board::HEADER), &scale).is_err());
        for args in [vec!["--velocity", "NaN"], vec!["--duration", "inf"], vec!["--modes", "0"],
            vec!["--scale", "input.csv", "--render", "input.csv"], vec!["--render"],
            vec!["--note", "69", "--note", "70"], vec!["--substeps", "17"]] {
            assert!(options(&args).is_err(), "accepted {args:?}");
        }
    }
    #[test]
    fn geometric_board_and_performance_options_compose_without_replacing_existing_controls() {
        let o = options(&["--scale", "strings.csv", "--board-geometry", "panel.fsb",
            "--board-band-hz", "300", "--performance", "score.csv", "--render", "piano.wav",
            "--sample-rate", "44100", "--duration", "2", "--modes", "40"]).unwrap();
        assert_eq!(o.board_geometry.as_deref(), Some("panel.fsb"));
        assert_eq!(o.performance.as_deref(), Some("score.csv"));
        assert_eq!(o.board_band_hz, 300.0);
        assert_eq!((o.sample_rate, o.duration, o.modes), (44_100, 2.0, 40));
        for args in [vec!["--board", "a.csv", "--board-geometry", "b.fsb"],
            vec!["--board-band-hz", "300"], vec!["--performance", "score.csv"],
            vec!["--render", "a.wav", "--performance", "score.csv", "--note", "69"],
            vec!["--render", "a.wav", "--performance", "score.csv", "--velocity", "1"],
            vec!["--render", "score.csv", "--performance", "score.csv"],
            vec!["--observer-gain", "NaN"],
            vec!["--board-geometry", "a.fsb", "--board-band-hz", "21600"]] {
            assert!(options(&args).is_err(), "accepted {args:?}");
        }
    }
    #[test]
    fn board_export_preserves_missing_measurements_and_reimports_the_admitted_subset() {
        let scale = vec![geometry::demonstration_scale().unwrap()[48]];
        let modes = board::demonstration();
        let text = write_board_for_scale(&modes, &scale);
        assert_eq!(text.lines().count(), modes.len() + 1);
        let roundtrip = board::read(&text, &[69]).unwrap();
        for (a, b) in roundtrip.iter().zip(&modes) {
            assert_eq!(a.frequency_hz, b.frequency_hz);
            assert_eq!(a.bridge[48], b.bridge[48]);
            assert_eq!(a.volume, b.volume);
        }
        assert!(board::read(&text, &[60, 69]).is_err());
    }
    #[test]
    fn source_preset_renders_and_exports_without_conflicting_board_inputs() {
        let o = options(&["--preset", "steinway-d", "--render", "d.wav", "--board-band-hz", "350",
            "--dump-geometry", "d.fsb", "--dump-obj", "d.obj", "--mesh-divisions", "12"]).unwrap();
        assert_eq!(o.preset.as_deref(), Some("steinway-d"));
        assert_eq!(o.mesh_divisions, 12);
        for args in [vec!["--preset", "unknown"],vec!["--preset", "steinway-d", "--board", "b.csv"],
            vec!["--preset", "steinway-d", "--board-geometry", "b.fsb", "--dump-obj", "d.obj"],vec!["--dump-obj", "d.obj"],
            vec!["--preset", "steinway-d", "--dump-obj", "d.obj", "--render", "d.obj"],
            vec!["--preset", "steinway-d", "--mesh-divisions", "3"]] {assert!(options(&args).is_err());}
    }
    #[test]
    fn physical_microphone_and_diagnostic_observer_are_not_silently_mixed() {
        let o=options(&["--preset","steinway-d","--render","d.wav","--microphone","0.5,1.2,0.8"]).unwrap();
        assert_eq!(o.microphone,Some([0.5,1.2,0.8]));assert!(!o.diagnostic_volume);
        let o=options(&["--preset","steinway-d","--render","d.wav","--diagnostic-volume","--observer-gain","500"]).unwrap();
        assert!(o.diagnostic_volume);
        for args in [vec!["--preset","steinway-d","--observer-gain","1000"],
            vec!["--render","d.wav","--microphone","0,0,1"],
            vec!["--preset","steinway-d","--render","d.wav","--microphone","0,0,NaN"],
            vec!["--preset","steinway-d","--render","d.wav","--microphone","0,0,-1"],
            vec!["--preset","steinway-d","--render","d.wav","--microphone","0,1"],
            vec!["--preset","steinway-d","--render","d.wav","--microphone","0,0,1","--diagnostic-volume"]] {
            assert!(options(&args).is_err());
        }
    }
    #[test]
    fn preset_tunes_tension_not_geometry_and_raw_source_is_available() {
        let o=options(&["--preset","steinway-d"]).unwrap();
        let tuned=selected_scale(None,&o).unwrap();let source=steinway_scale::courses().unwrap();
        for (a,b) in tuned.iter().zip(&source) {
            assert_eq!(a.length_m,b.length_m);assert_eq!(a.linear_density_kg_m,b.linear_density_kg_m);
            assert_eq!(a.flexural_rigidity_nm2,b.flexural_rigidity_nm2);
            let target=440.0*2.0f64.powf((f64::from(a.midi)-69.0)/12.0);
            assert!((a.partial_hz(1,a.tension_n)/target-1.0).abs()<1e-10);
        }
        let raw=options(&["--preset","steinway-d","--raw-tensions"]).unwrap();
        assert_eq!(selected_scale(None,&raw).unwrap(),source);
        let imported=options(&["--preset","steinway-d","--scale","measured.csv"]).unwrap();
        assert_eq!(selected_scale(Some(&geometry::write_scale(&source)),&imported).unwrap(),source);
        for args in [vec!["--concert-pitch","NaN"],vec!["--concert-pitch","400"],
            vec!["--raw-tensions"],vec!["--preset","steinway-d","--raw-tensions","--concert-pitch","442"]] {
            assert!(options(&args).is_err());
        }
    }
    #[test]
    fn source_hammer_cards_are_used_by_the_render_preparation_path() {
        let mut o=options(&["--preset","steinway-d"]).unwrap();o.modes=12;
        let scale=selected_scale(None,&o).unwrap();
        let mut piano=prepare_instrument(vec![scale[48]],&board::demonstration(),&o).unwrap();
        piano.note_on(69,2.0).unwrap();
        for _ in 0..1500 {assert!(piano.step().unwrap().is_finite());}
        assert!(piano.accounting.felt_loss_j>0.0);
        assert!(piano.accounting.felt_relaxation_loss_j>0.0);
        assert!((piano.accounting.input_work_j-piano.energy_j()-piano.accounting.dissipated_j()).abs()<1e-7);
    }
    #[test]
    fn jack_performance_reaches_preset_mechanics_without_velocity_substitution() {
        let mut o=options(&["--preset","steinway-d"]).unwrap();o.modes=12;
        let c=selected_scale(None,&o).unwrap()[48];
        let mut piano=prepare_instrument(vec![c],&board::demonstration(),&o).unwrap();
        let mut score=performance::Performance::read("sample,event,key,value\n0,jack_legato,69,30\n",&[69],2400).unwrap();
        for sample in 0..2400 {score.dispatch(sample,&mut piano).unwrap();piano.step().unwrap();}
        assert!(piano.accounting.shank_loss_j>0.0);assert!(piano.accounting.felt_loss_j>0.0);
        assert!((piano.accounting.input_work_j-piano.energy_j()-piano.accounting.dissipated_j()).abs()<1e-7);
    }
    #[test]
    fn supplied_hammers_compose_with_presets_and_protect_the_input_file() {
        let o = options(&["--preset", "steinway-d", "--hammers", "felt.fsh", "--render", "d.wav"]).unwrap();
        assert_eq!(o.hammers.as_deref(), Some("felt.fsh"));
        assert!(options(&["--hammers", "felt.fsh"]).is_err());
        assert!(options(&["--hammers", "felt.fsh", "--render", "felt.fsh"]).is_err());
        let course = selected_scale(None, &o).unwrap()[48];
        let text = "frankensim-hammer-materials-v1\nfelt,69,400000,0.2,2.5,3.2,0.25,0.8,2500000\n";
        let mut piano = prepare_instrument_with_hammers(vec![course], &board::demonstration(),
            &o, Some(text)).unwrap();
        // An imported elastic card must not resurrect the preset Prony tails,
        // while the preset's geometry-derived jack/shank remains available.
        piano.jack_on(69, 30.0, 0.1).unwrap();
        for _ in 0..1500 { assert!(piano.step().unwrap().is_finite()); }
        assert_eq!(piano.accounting.felt_relaxation_loss_j, 0.0);
        assert!(piano.accounting.shank_loss_j > 0.0);
        assert!(prepare_instrument_with_hammers(vec![course], &board::demonstration(),
            &o, Some("invalid material")).is_err());
    }
}
