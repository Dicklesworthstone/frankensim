//! Geometry-first grand-piano composition of existing FrankenSim owners.
//! `cargo run --release -p fs-couple --example grand_piano -- --render piano.wav`
//! Defaults are ESTIMATED, not a measured Steinway digital twin.
mod geometry;
mod linear;
mod board;
mod felt;
mod engine;

const USAGE: &str = "grand_piano [--render piano.wav] [--scale strings.csv] [--board board.csv]
    [--note 21..108] [--velocity m/s] [--duration seconds]
    [--sample-rate Hz] [--substeps 1..16] [--modes 1..256]
    [--dump-scale strings.csv] [--dump-board board.csv]
Imports are used by the renderer, not just inspected. Input CSV values are not
independently certified as Steinway measurements. Omitted inputs are estimates.
--note performs a single-key study; otherwise the demo also plays a chord of
available keys. --velocity overrides the three demo hammer launch speeds.
Velocity is POST-ESCAPEMENT hammer velocity, not MIDI velocity or key motion.";

#[derive(Debug)]
struct Options {
    render: Option<String>, scale: Option<String>, board: Option<String>,
    dump_scale: Option<String>, dump_board: Option<String>,
    note: Option<u8>, velocity: Option<f64>, duration: f64,
    sample_rate: u32, substeps: usize, modes: usize, help: bool,
}
impl Default for Options {
    fn default() -> Self {
        Self { render: None, scale: None, board: None, dump_scale: None,
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
            let value = args.next().ok_or_else(|| format!("missing value for {flag}"))?;
            let invalid = || format!("invalid value for {flag}: {value}");
            match flag.as_str() {
                "--render" => options.render = Some(value.clone()),
                "--scale" => options.scale = Some(value.clone()),
                "--board" => options.board = Some(value.clone()),
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
            || !(1..=16).contains(&options.substeps) || !(1..=256).contains(&options.modes) {
            return Err("render control outside its finite admitted range".into());
        }
        // Do not overwrite the very measurements that a render was asked to use.
        let inputs = [options.scale.as_ref(), options.board.as_ref()];
        let outputs = [options.render.as_ref(), options.dump_scale.as_ref(), options.dump_board.as_ref()];
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
}

fn load_scale(text: Option<&str>) -> Result<Vec<geometry::Course>, String> {
    match text { Some(text) => geometry::read_scale(text), None => geometry::demonstration_scale() }
}
fn load_board(text: Option<&str>, scale: &[geometry::Course]) -> Result<Vec<linear::BoardMode>, String> {
    match text {
        Some(text) => board::read(text, &scale.iter().map(|c| c.midi).collect::<Vec<_>>()),
        None => Ok(board::demonstration()),
    }
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
    options: &Options) -> Result<(), String> {
    let key = study_key(&scale, options.note)?;
    let chord: Vec<u8> = if options.note.is_some() { Vec::new() } else {
        [48, 60, 64, 67].into_iter().filter(|key| scale.iter().any(|c| c.midi == *key)).collect()
    };
    let mut piano = engine::Instrument::new(scale, modes, options.sample_rate,
        options.substeps, options.modes, true)?;
    piano.set_sustain(1.0).map_err(|e| e.to_string())?;
    let rate = piano.sample_rate();
    let count = (options.duration * f64::from(rate)).round() as u32;
    let mut pressure = Vec::with_capacity(count as usize);
    let start = std::time::Instant::now();
    let mut peak: f64 = 0.0;
    for sample in 0..count {
        for (second, velocity) in [(0, 0.6), (1, 2.0), (2, 4.5)] {
            if sample == rate * second {
                piano.note_on(key, options.velocity.unwrap_or(velocity)).map_err(|e| e.to_string())?;
            }
        }
        if [rate / 2, rate + rate / 2, rate * 2 + rate / 2].contains(&sample) {
            piano.note_off(key).map_err(|e| e.to_string())?;
        }
        if sample == rate * 3 {
            for &key in &chord { piano.note_on(key, options.velocity.unwrap_or(2.5)).map_err(|e| e.to_string())?; }
        }
        if sample == rate * 4 {
            for &key in &chord { piano.note_off(key).map_err(|e| e.to_string())?; }
        }
        if sample == rate * 4 + rate / 2 { piano.set_sustain(0.5).map_err(|e| e.to_string())?; }
        if sample == rate * 5 { piano.set_sustain(0.0).map_err(|e| e.to_string())?; }
        // Explicit diagnostic observer gain [Pa / (m^3/s)]. This is NOT a
        // measured radiation transfer or a claim of calibrated acoustic SPL.
        let p = 10_000.0 * piano.step().map_err(|e| format!("sample {sample}: {e}"))?;
        peak = peak.max(p.abs());
        pressure.push(p);
    }
    let elapsed = start.elapsed().as_secs_f64();
    let seconds = f64::from(count) / f64::from(rate);
    let (wav, clips) = fs_couple::pcm_wav::encode_pcm16_wav(&pressure, rate, 2.0).map_err(|e| e.to_string())?;
    std::fs::write(path, wav).map_err(|e| format!("{path}: {e}"))?;
    println!("Diagnostic volume-velocity observer, gain 10000 Pa/(m^3/s); no peak normalization.");
    println!("{} string modes; {} board modes; {} above-band duplex segments omitted from dynamic retention (static attachment retained).",
        piano.bank.modes.len(), piano.bank.board_count, piano.bank.omitted_duplex_modes);
    println!("{seconds:.6} s audio rendered in {elapsed:.6} s; wall/audio ratio {:.4}; peak {peak:.6} Pa-equivalent; {clips} PCM clips.", elapsed / seconds);
    println!("Input {:.9} J; stored {:.9} J; component losses {:.9} J; closure {:.3e} J; worst substep defect {:.3e} J.",
        piano.accounting.input_work_j, piano.energy_j(), piano.accounting.dissipated_j(),
        piano.accounting.input_work_j - piano.energy_j() - piano.accounting.dissipated_j(), piano.accounting.max_balance_error_j);
    Ok(())
}

fn run() -> Result<(), String> {
    let options = Options::parse(&std::env::args().skip(1).collect::<Vec<_>>())?;
    if options.help { println!("{USAGE}"); return Ok(()); }
    let read = |path: &String| std::fs::read_to_string(path).map_err(|e| format!("{path}: {e}"));
    let scale_text = options.scale.as_ref().map(read).transpose()?;
    let board_text = options.board.as_ref().map(read).transpose()?;
    let scale = load_scale(scale_text.as_deref())?;
    let modes = load_board(board_text.as_deref(), &scale)?;
    study_key(&scale, options.note)?;
    println!("String scale: {}.", options.scale.as_deref().unwrap_or("ESTIMATED demonstration"));
    println!("Soundboard: {}.", options.board.as_deref().unwrap_or("AUTHORED illustrative modes"));
    println!("Source authority belongs to the inputs, not the model name; imported files are not independently certified measurements.");
    if let Some(path) = &options.dump_scale {
        let source = options.scale.as_deref().unwrap_or("ESTIMATED demonstration; not measured Steinway geometry");
        std::fs::write(path, format!("# Source: {source}\n{}", geometry::write_scale(&scale))).map_err(|e| e.to_string())?;
    }
    if let Some(path) = &options.dump_board {
        let source = options.board.as_deref().unwrap_or("AUTHORED illustrative board; not measured Steinway geometry");
        std::fs::write(path, format!("# Source: {source}\n{}", board::write(&modes))).map_err(|e| e.to_string())?;
    }
    if let Some(path) = &options.render { return render(path, scale, &modes, &options); }
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
mod render_tests {
    use super::*;
    fn options(args: &[&str]) -> Result<Options, String> {
        Options::parse(&args.iter().map(|s| (*s).to_owned()).collect::<Vec<_>>())
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
}
